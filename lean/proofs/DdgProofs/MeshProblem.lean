import DdgProofs.Certificate
import DdgProofs.Problem
import DDG.Mesh.GaussBonnet
import DDG.Mesh.TriMesh
import DDG.Laplace.Cotan

/-!
# The canonical problem of a mesh, and vectors on all vertices

`ConeProblem ι` is stated on the unknowns `ι` only. The Rust code (and `CLAUDE.md`) work with
one value per mesh vertex, zero at the vertices that are not unknowns (the boundary). This file
ties the two together:

* `TriMesh.coneProblem M p wE lam site` is the canonical problem of the mesh `M` whose unknowns
  are the vertices in `p` (the interior vertices; all vertices on a closed mesh): `L` is the
  block of the cotan Laplacian, `Ω` the angle defects, `m` the lumped masses divided by the
  total area of the mesh.
* `coneProblem_stateResidual`, `coneProblem_adjoint`: for vectors on all vertices that vanish off
  `p`, the residuals of `ConeProblem` are the full-mesh expressions `(L u)_v − μ_v + Ω_v` and
  `(L φ)_v` at the unknowns — the quantities the certificate code evaluates face by face.
* `coneProblem_feasible_iff`: feasibility is the constraint set of `CLAUDE.md`,
  `(L u)_v = μ_v − Ω_v` and `μ_v` admissible at the unknowns, with `u = 0` elsewhere.
* `coneProblem_objective`: the objective is `½ Σ_{v ∈ p} wE_v (m_v/A) u_v² + λ Σ_{v ∈ p} cost_v(μ_v)`.
* `coneProblem_transpose`, `coneProblem_posDef`, `coneProblem_mulVec_one`: the matrix is
  symmetric; positive definite if the mesh is connected and some vertex is not an unknown;
  annihilates constants if every vertex is an unknown.

So every theorem about `ConeProblem` (optimality, certificate, solver) applies to per-vertex
arrays `u, φ, μ` that are zero off the unknowns, a condition that can be tested exactly.
-/

noncomputable section

namespace DDG

open Matrix

namespace TriMesh

variable {V F : Type*} [DecidableEq V] [Fintype V] [Fintype F] (M : TriMesh V F)

/-- The canonical discrete problem of the mesh `M` on the unknowns `{v // p v}`: the block of
the cotan Laplacian, the angle defects, and the lumped masses divided by the total area. -/
def coneProblem (p : V → Prop) [DecidablePred p] (wE : V → ℝ) (lam : ℝ) (site : V → Site) :
    ConeProblem {v // p v} where
  L := M.laplacian.submatrix Subtype.val Subtype.val
  Ω := fun i => M.angleDefect i
  m := fun i => M.mass i / M.totalArea
  wE := fun i => wE i
  lam := lam
  site := fun i => site i

variable (p : V → Prop) [DecidablePred p] (wE : V → ℝ) (lam : ℝ) (site : V → Site)

/-- The matrix of the mesh problem is symmetric. -/
theorem coneProblem_transpose : (M.coneProblem p wE lam site).Lᵀ = (M.coneProblem p wE lam site).L :=
  (M.laplacian_isSymm.submatrix _ : _)

/-- **Dirichlet case**: on a connected valid mesh with a vertex that is not an unknown, the
matrix of the mesh problem is positive definite. -/
theorem coneProblem_posDef (hM : M.Valid) (hconn : M.Connected) (hB : ∃ v, ¬p v) :
    (M.coneProblem p wE lam site).L.PosDef :=
  M.laplacian_interior_posDef hM hconn p hB

/-- **The state residual is the full-mesh expression.** If `u` vanishes off the unknowns, the
state residual of the mesh problem at an unknown `v` is `(L u)_v − μ_v + Ω_v` with the full cotan
Laplacian applied to the per-vertex vector. -/
theorem coneProblem_stateResidual {u : V → ℝ} (hu : ∀ v, ¬p v → u v = 0) (μ : V → ℝ)
    (i : {v // p v}) :
    (M.coneProblem p wE lam site).stateResidual (fun j => u j) (fun j => μ j) i =
      (M.laplacian *ᵥ u) i - μ i + M.angleDefect i := by
  simp only [ConeProblem.stateResidual, coneProblem]
  rw [← M.laplacian_mulVec_of_zero_off p u hu i]

/-- **The adjoint product is the full-mesh expression**: if `φ` vanishes off the unknowns,
`(Lᵀφ)_v = (L φ)_v` at an unknown `v`, with the full cotan Laplacian. -/
theorem coneProblem_adjoint {φ : V → ℝ} (hφ : ∀ v, ¬p v → φ v = 0) (i : {v // p v}) :
    ((M.coneProblem p wE lam site).Lᵀ *ᵥ fun j : {v // p v} => φ j) i = (M.laplacian *ᵥ φ) i := by
  rw [coneProblem_transpose]
  simp only [coneProblem]
  rw [← M.laplacian_mulVec_of_zero_off p φ hφ i]

/-- **Feasibility is the constraint set of `CLAUDE.md`.** For per-vertex vectors with `u = 0`
off the unknowns: the restriction is feasible for the mesh problem iff `(L u)_v = μ_v − Ω_v` and
`μ_v` is admissible at every unknown `v`. -/
theorem coneProblem_feasible_iff {u : V → ℝ} (hu : ∀ v, ¬p v → u v = 0) (μ : V → ℝ) :
    (M.coneProblem p wE lam site).Feasible (fun j => u j) (fun j => μ j) ↔
      (∀ v, p v → (M.laplacian *ᵥ u) v = μ v - M.angleDefect v) ∧
        ∀ v, p v → μ v ∈ (site v).feasible := by
  unfold ConeProblem.Feasible
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨fun v hv => ?_, fun v hv => h2 ⟨v, hv⟩⟩
    have := congrFun h1 ⟨v, hv⟩
    simp only [coneProblem, Pi.sub_apply] at this
    rw [M.laplacian_mulVec_of_zero_off p u hu ⟨v, hv⟩]
    exact this
  · rintro ⟨h1, h2⟩
    refine ⟨?_, fun i => h2 i i.2⟩
    funext i
    simp only [coneProblem, Pi.sub_apply]
    rw [← M.laplacian_mulVec_of_zero_off p u hu i]
    exact h1 i i.2

/-- **The objective in per-vertex terms**:
`½ Σ_{v ∈ p} wE_v (m_v/A) u_v² + λ Σ_{v ∈ p} cost_v(μ_v)`, `A` the total area of the mesh. -/
theorem coneProblem_objective (u μ : V → ℝ) :
    (M.coneProblem p wE lam site).objective (fun j => u j) (fun j => μ j) =
      1 / 2 * ∑ v ∈ Finset.univ.filter p, wE v * (M.mass v / M.totalArea) * u v ^ 2 +
        lam * ∑ v ∈ Finset.univ.filter p, (site v).cost (μ v) := by
  simp only [ConeProblem.objective, ConeProblem.distortion, ConeProblem.coneCost, coneProblem]
  rw [← Finset.sum_subtype (Finset.univ.filter p) (by simp)
      (fun v => wE v * (M.mass v / M.totalArea) * u v ^ 2),
    ← Finset.sum_subtype (Finset.univ.filter p) (by simp) (fun v => (site v).cost (μ v))]

/-- **Closed case**: if every vertex is an unknown, the matrix of the mesh problem annihilates
constants (and so has zero column sums, being symmetric). -/
theorem coneProblem_mulVec_one (hp : ∀ v, p v) :
    (M.coneProblem p wE lam site).L *ᵥ (fun _ => (1 : ℝ)) = 0 := by
  funext i
  have h := M.laplacian_mulVec_of_zero_off p (fun _ => (1 : ℝ)) (fun v hv => absurd (hp v) hv) i
  simp only [coneProblem]
  rw [← h, M.laplacian_mulVec_one]
  rfl

/-- The mesh problem is well-formed when the weights are: `wE ≥ 0`, `λ ≥ 0`, well-formed sites,
on a valid mesh (which makes the masses nonnegative). -/
theorem coneProblem_wellFormed (hM : M.Valid) (hwE : ∀ v, p v → 0 ≤ wE v) (hlam : 0 ≤ lam)
    (hsite : ∀ v, p v → (site v).WellFormed) : (M.coneProblem p wE lam site).WellFormed where
  wE_nonneg := fun i => hwE i i.2
  m_nonneg := fun i => by
    simp only [coneProblem]
    have hA : 0 ≤ M.totalArea :=
      Finset.sum_nonneg fun f _ => (triArea_pos (hM.tri f)).le
    have hm : 0 ≤ M.mass i := by
      unfold mass vertexArea
      refine div_nonneg (Finset.sum_nonneg fun f _ => ?_) (by norm_num)
      split_ifs
      · exact (triArea_pos (hM.tri f)).le
      · exact le_rfl
    exact div_nonneg hm hA
  lam_nonneg := hlam
  site := fun i => hsite i i.2

end TriMesh

end DDG
