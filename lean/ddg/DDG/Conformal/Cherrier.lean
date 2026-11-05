import DDG.Conformal.Yamabe
import DDG.Laplace.Poisson
import DDG.Conformal.BoundaryCurve
import DDG.Laplace.Cotan

/-!
# Boundary First Flattening: the linearized curvatures, Cherrier's formula, the Hilbert transform

`crates/ddg-bff` (Sawhney & Crane, *Boundary First Flattening*, ACM TOG 2017; "BFF"), on top of
the Poisson problems of `DdgProofs.Laplace.Poisson`.

* **`bff-exterior-angle-sum`** (BFF Proposition 1): `gauss_bonnet_conformal`: for every `u`,
  `Σ_I (Ω + L u) + Σ_B (k + L u) = 2πχ`. Hence the exterior angles `k̃ = k + (L u)_B` of the
  target curve sum to `2πχ` when `L u = −Ω` in the interior (`sum_exteriorAngle`), and to
  `2πχ − Σ_I μ` with cones (`sum_exteriorAngle_cones`).
* **`bff-cherrier-boundary-curvature`** (BFF Eq. 12): `hasDerivAt_geodesicCurvature_conformalScale`:
  under `ℓ̃_ij = e^{t(v_i+v_j)/2} ℓ_ij` the boundary curvature `k_i = π − Σθ` has derivative
  `(L v)_i` at `t = 0`. So `k̃ = k + (L u)_B` is the first-order change of the boundary curvature;
  the code *defines* the target exterior angles by it.

* **`bff-hilbert-transform`** (BFF Eq. 17): `hilbertTransform_cotan` puts the statements of
  `DDG/Conformal/BoundaryCurve.lean` on a mesh: the Neumann problem for the harmonic conjugate is
  solvable and its solutions are the minimizers of the conformal energy.

Signs: `L` is positive semidefinite and `u` is the log scale factor (`CLAUDE.md`); BFF's `A` is
this `L` and its Eq. 12–14 read `A u = −(K − C)` in the interior and `h = −(A u)_B` on the
boundary for the source `φ = −Ω`.
-/

noncomputable section

namespace DDG

open Matrix Real

namespace TriMesh

variable {V F : Type*} [DecidableEq V] [Fintype V] [Fintype F] (M : TriMesh V F)

/-! ## Gauss–Bonnet for the linearized curvatures; the exterior angles sum to `2πχ` -/

/-- **Gauss–Bonnet for the first-order curvatures of a conformally scaled metric.** For a valid
manifold mesh and *every* `u`:
`Σ_{v interior} (Ω_v + (L u)_v) + Σ_{v boundary} (k_v + (L u)_v) = 2πχ`.
(The columns of `L` sum to zero, so the correction terms cancel in total.)

[formula: bff-exterior-angle-sum] -/
theorem gauss_bonnet_conformal (hM : M.Valid) (hman : M.Manifold) (u : V → ℝ) :
    ∑ v ∈ M.boundaryVerticesᶜ, (M.angleDefect v + (M.laplacian *ᵥ u) v) +
        ∑ v ∈ M.boundaryVertices, (M.geodesicCurvature v + (M.laplacian *ᵥ u) v) =
      2 * π * (M.eulerChar : ℝ) := by
  have h0 : ∑ v, (M.laplacian *ᵥ u) v = 0 :=
    Solver.sum_eq_zero_of_mulVec_eq M.laplacian_vecMul_one rfl
  have hsplit : ∑ v ∈ M.boundaryVerticesᶜ, (M.laplacian *ᵥ u) v +
      ∑ v ∈ M.boundaryVertices, (M.laplacian *ᵥ u) v = ∑ v, (M.laplacian *ᵥ u) v := by
    rw [add_comm]; exact Finset.sum_add_sum_compl _ _
  have hgb := M.gauss_bonnet hM hman
  simp only [Finset.sum_add_distrib]
  linarith

/-- **The exterior angles of the target curve sum to `2πχ`** (BFF Proposition 1; `2π` for a
disk, `χ = 1`). Let `u` satisfy `(L u)_v = −Ω_v` at the interior vertices (the scale factors of
a flat target metric, any boundary values). Then the exterior angles
`k̃_v = k_v + (L u)_v` (`= k_v − h_v` with the Neumann data `h = −(L u)_B` of the source `−Ω`)
satisfy `Σ_{v boundary} k̃_v = 2πχ`.

[formula: bff-exterior-angle-sum] -/
theorem sum_exteriorAngle (hM : M.Valid) (hman : M.Manifold) {u : V → ℝ}
    (hu : ∀ v ∈ M.boundaryVerticesᶜ, (M.laplacian *ᵥ u) v = -M.angleDefect v) :
    ∑ v ∈ M.boundaryVertices, (M.geodesicCurvature v + (M.laplacian *ᵥ u) v) =
      2 * π * (M.eulerChar : ℝ) := by
  have h := M.gauss_bonnet_conformal hM hman u
  have h0 : ∑ v ∈ M.boundaryVerticesᶜ, (M.angleDefect v + (M.laplacian *ᵥ u) v) = 0 :=
    Finset.sum_eq_zero fun v hv => by rw [hu v hv]; ring
  linarith

/-- **With cones**: if `(L u)_v = μ_v − Ω_v` at the interior vertices, the exterior angles sum to
`2πχ − Σ_{v interior} μ_v`.

[formula: bff-exterior-angle-sum] -/
theorem sum_exteriorAngle_cones (hM : M.Valid) (hman : M.Manifold) {u μ : V → ℝ}
    (hu : ∀ v ∈ M.boundaryVerticesᶜ, (M.laplacian *ᵥ u) v = μ v - M.angleDefect v) :
    ∑ v ∈ M.boundaryVertices, (M.geodesicCurvature v + (M.laplacian *ᵥ u) v) =
      2 * π * (M.eulerChar : ℝ) - ∑ v ∈ M.boundaryVerticesᶜ, μ v := by
  have h := M.gauss_bonnet_conformal hM hman u
  have h0 : ∑ v ∈ M.boundaryVerticesᶜ, (M.angleDefect v + (M.laplacian *ᵥ u) v) =
      ∑ v ∈ M.boundaryVerticesᶜ, μ v :=
    Finset.sum_congr rfl fun v hv => by rw [hu v hv]; ring
  linarith

/-! ## The discrete Cherrier formula -/

/-- **First-order change of the boundary curvature under a conformal change of metric.** Under
`ℓ̃_ij = e^{t(v_i + v_j)/2}·ℓ_ij`, the curvature `k_i = π − Σθ` at a vertex `i` has derivative
`(L v)_i` at `t = 0` — at a boundary vertex this is minus the Neumann data of `v` for the cotan
Laplacian. BFF Eq. 12 (`k̃ = k − h`) is this first-order formula; in the smooth setting the
corresponding identity (Cherrier) is exact.

[formula: bff-cherrier-boundary-curvature] -/
theorem hasDerivAt_geodesicCurvature_conformalScale (hM : M.Valid) (v : V → ℝ) (i : V) :
    HasDerivAt (fun t : ℝ => (M.conformalScale (fun j => t * v j)).geodesicCurvature i)
      ((M.laplacian *ᵥ v) i) 0 := by
  have h := (M.hasDerivAt_angleSum_conformalScale hM v i).const_sub π
  simpa [geodesicCurvature] using h

/-! ## The harmonic conjugate (Hilbert transform) on a mesh -/

/-- **The Neumann problem of the Hilbert transform is solvable, and its solutions are the
minimizers of the conformal energy.** On a connected valid mesh, for a boundary loop `β` with
successor `next`, an injective map `emb` from boundary positions to vertices, and any `a`:
there is a `b` with `L b = h`, where `h_j = ½(a_{j−1} − a_{j+1})` at the boundary vertex `emb j`
and `0` at every other vertex (`hilbertVector`, the vector `hilbert_transform` returns); and `b`
solves this system iff it minimizes `b' ↦ ½(aᵀLa + b'ᵀLb') − 𝒜(a, b')`, `𝒜` the signed area of
the image of the boundary polygon (`loopConformalEnergy`). Two solutions differ by a constant
(`sub_const_of_laplacian_mulVec_eq`).

[formula: bff-hilbert-transform] -/
theorem hilbertTransform_cotan {β : Type*} [Fintype β] (hM : M.Valid) (hconn : M.Connected)
    (next : Equiv.Perm β) {emb : β → V} (hemb : Function.Injective emb) (a : V → ℝ) :
    (∃ b, M.laplacian *ᵥ b = hilbertVector next emb a) ∧
      ∀ b, (∀ b', loopConformalEnergy next emb M.laplacian a b ≤
            loopConformalEnergy next emb M.laplacian a b') ↔
        M.laplacian *ᵥ b = hilbertVector next emb a := by
  refine ⟨?_, fun b => hilbertVector_minimizes_conformalEnergy_iff next hemb
    (M.laplacian_posSemidef hM) a b⟩
  rw [M.exists_laplacian_mulVec_eq_iff hM hconn]
  exact sum_hilbertRhs_vertices hemb (fun j => hilbertVector_emb next hemb a j)
    (fun v hv => hilbertVector_of_notMem next emb a hv)

end TriMesh

end DDG
