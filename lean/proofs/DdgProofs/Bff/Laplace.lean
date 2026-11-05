import DdgProofs.Yamabe
import DdgProofs.Solver
import DdgProofs.Bff.Curve

/-!
# The Poisson problems of Boundary First Flattening

`crates/ddg-bff` solves three linear problems with the cotan Laplacian `L` of a mesh (Sawhney &
Crane, *Boundary First Flattening*, ACM TOG 2017; "BFF"). With `F` the set of free vertices and
`X` its complement (the vertices whose value is prescribed):

* **`bff-harmonic-extension`** (`PoissonSystem::solve`): `L_FF a_F = φ_F − L_FX a_X`.
  `existsUnique_poisson`: on a connected valid mesh with `X ≠ ∅` there is exactly one `a` with
  the prescribed values on `X` and `(L a)_v = φ_v` on `F`. With `φ = 0` this is the harmonic
  extension of boundary values (BFF Sec. 4.6, Alg. 1 line 11).
* **`bff-dirichlet-to-neumann`** (`dirichlet_to_neumann`, BFF Eq. 14, Alg. 3): from the solution
  `a` of the Dirichlet problem, the Neumann data `h_B = φ_B − (L a)_B`. `neumannData_spec`:
  `L a = φ − [0; h]`, i.e. `a` solves the Neumann problem with this `h`;
  `sum_neumannData`: `Σ_B h = Σ_V φ`, the compatibility condition of that Neumann problem.
* **`bff-neumann-to-dirichlet`** (`neumann_to_dirichlet`, BFF Eq. 5, Alg. 4): `L a = b`.
  `sum_eq_zero_of_laplacian_mulVec_eq`: solvable only if `Σ b = 0`; `existsUnique_pinned`:
  if `Σ b = 0` there is exactly one solution with `a(v₀) = 0` (the code pins a vertex and drops
  its equation); `existsUnique_neumann`: exactly one solution with `Σ m a = 0` (the code then
  subtracts the mass-weighted mean); any two solutions differ by a constant.
* **`bff-exterior-angle-sum`** (BFF Proposition 1): `gauss_bonnet_conformal`: for every `u`,
  `Σ_I (Ω + L u) + Σ_B (k + L u) = 2πχ`. Hence the exterior angles `k̃ = k + (L u)_B` of the
  target curve sum to `2πχ` when `L u = −Ω` in the interior (`sum_exteriorAngle`), and to
  `2πχ − Σ_I μ` with cones (`sum_exteriorAngle_cones`).
* **`bff-cherrier-boundary-curvature`** (BFF Eq. 12): `hasDerivAt_geodesicCurvature_conformalScale`:
  under `ℓ̃_ij = e^{t(v_i+v_j)/2} ℓ_ij` the boundary curvature `k_i = π − Σθ` has derivative
  `(L v)_i` at `t = 0`. So `k̃ = k + (L u)_B` is the first-order change of the boundary curvature;
  the code *defines* the target exterior angles by it.

* **`bff-hilbert-transform`** (BFF Eq. 17): `hilbertTransform_cotan` puts the statements of
  `DdgProofs/Bff/Curve.lean` on a mesh: the Neumann problem for the harmonic conjugate is
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

/-! ## The Dirichlet (Poisson) problem -/

/-- **Existence and uniqueness for the discrete Poisson problem with prescribed values.** On a
connected valid mesh, let `p` be the set of free vertices and let at least one vertex be
prescribed. For every source `φ` and all prescribed values `g` there is exactly one `a` with
`a = g` off `p` and `(L a)_v = φ_v` on `p`.

This is the system `L_FF a_F = φ_F − L_FX a_X` that `PoissonSystem::solve` factorizes; the
interior block is positive definite (`laplacian_interior_posDef`).

[formula: bff-harmonic-extension] -/
theorem existsUnique_poisson (hM : M.Valid) (hconn : M.Connected) (p : V → Prop)
    [DecidablePred p] (hB : ∃ v, ¬p v) (φ g : V → ℝ) :
    ∃! a : V → ℝ, (∀ v, ¬p v → a v = g v) ∧ ∀ v, p v → (M.laplacian *ᵥ a) v = φ v := by
  set A := M.laplacian.submatrix (Subtype.val : {v // p v} → V) Subtype.val with hA
  have hpd : A.PosDef := M.laplacian_interior_posDef hM hconn p hB
  have hdet : IsUnit A.det := (Matrix.isUnit_iff_isUnit_det _).mp hpd.isUnit
  set g₀ : V → ℝ := fun v => if p v then 0 else g v with hg₀
  set r : {v // p v} → ℝ := fun i => φ i - (M.laplacian *ᵥ g₀) i with hr
  set x : {v // p v} → ℝ := A⁻¹ *ᵥ r with hx
  set z : V → ℝ := fun v => if h : p v then x ⟨v, h⟩ else 0 with hz
  have hzoff : ∀ v, ¬p v → z v = 0 := fun v hv => by simp [hz, hv]
  have hzx : (fun j : {v // p v} => z j) = x := by funext j; simp [hz, j.2]
  refine ⟨g₀ + z, ⟨fun v hv => ?_, fun v hv => ?_⟩, ?_⟩
  · simp [hg₀, hv, hzoff v hv]
  · have h1 := M.laplacian_mulVec_of_zero_off p z hzoff ⟨v, hv⟩
    rw [hzx] at h1
    have h2 : (A *ᵥ x) ⟨v, hv⟩ = r ⟨v, hv⟩ := by
      rw [hx, Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hdet, Matrix.one_mulVec]
    rw [Matrix.mulVec_add]
    simp only [Pi.add_apply]
    rw [h1, h2]
    simp [hr]
  · rintro b ⟨hb1, hb2⟩
    set d : V → ℝ := b - (g₀ + z) with hd
    have hdoff : ∀ v, ¬p v → d v = 0 := fun v hv => by
      simp [hd, hb1 v hv, hg₀, hv, hzoff v hv]
    have hLd : ∀ v, p v → (M.laplacian *ᵥ d) v = 0 := by
      intro v hv
      have h1 := M.laplacian_mulVec_of_zero_off p z hzoff ⟨v, hv⟩
      rw [hzx] at h1
      have h2 : (A *ᵥ x) ⟨v, hv⟩ = r ⟨v, hv⟩ := by
        rw [hx, Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hdet, Matrix.one_mulVec]
      rw [hd, Matrix.mulVec_sub, Matrix.mulVec_add]
      simp only [Pi.sub_apply, Pi.add_apply]
      rw [hb2 v hv, h1, h2]
      simp [hr]
    have hAd : A *ᵥ (fun j : {v // p v} => d j) = 0 := by
      funext i
      rw [← M.laplacian_mulVec_of_zero_off p d hdoff i]
      exact hLd i i.2
    have hd0 : (fun j : {v // p v} => d j) = 0 := by
      have := congrArg (fun y => A⁻¹ *ᵥ y) hAd
      simpa [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hdet] using this
    funext v
    by_cases hv : p v
    · have := congrFun hd0 ⟨v, hv⟩
      simp only [hd, Pi.sub_apply, Pi.zero_apply] at this
      linarith
    · have := hdoff v hv
      simp only [hd, Pi.sub_apply] at this
      linarith

/-! ## Dirichlet to Neumann -/

/-- The Neumann data of a function `a` for the source `φ`: `h_v = φ_v − (L a)_v` at the
prescribed (boundary) vertices, `0` at the free ones (`dirichlet_to_neumann` returns the boundary
part). -/
def neumannData (p : V → Prop) [DecidablePred p] (φ a : V → ℝ) : V → ℝ :=
  fun v => if p v then 0 else φ v - (M.laplacian *ᵥ a) v

/-- **The partitioned system.** If `(L a)_v = φ_v` at the free vertices, then with the Neumann
data `h` of `a`: `L a = φ − h` at *every* vertex. So the solution of the Dirichlet problem is a
solution of the Neumann problem with the data it produces.

[formula: bff-dirichlet-to-neumann] -/
theorem neumannData_spec (p : V → Prop) [DecidablePred p] (φ a : V → ℝ)
    (h : ∀ v, p v → (M.laplacian *ᵥ a) v = φ v) :
    M.laplacian *ᵥ a = fun v => φ v - M.neumannData p φ a v := by
  funext v
  by_cases hv : p v
  · simp [neumannData, hv, h v hv]
  · simp [neumannData, hv]

/-- **The Neumann data are compatible**: if `(L a)_v = φ_v` at the free vertices then
`Σ_{v prescribed} h_v = Σ_{v} φ_v` (because the columns of `L` sum to zero). This is the condition
under which the Neumann problem `L a' = φ − h` can be solved.

[formula: bff-dirichlet-to-neumann] -/
theorem sum_neumannData (p : V → Prop) [DecidablePred p] (φ a : V → ℝ)
    (h : ∀ v, p v → (M.laplacian *ᵥ a) v = φ v) :
    ∑ v, M.neumannData p φ a v = ∑ v, φ v := by
  have h1 := Solver.sum_eq_zero_of_mulVec_eq M.laplacian_vecMul_one
    (M.neumannData_spec p φ a h)
  rw [Finset.sum_sub_distrib] at h1
  linarith

/-! ## The Neumann problem -/

/-- A system `L a = b` with the cotan Laplacian is solvable only if `Σ b = 0`.

[formula: bff-neumann-to-dirichlet] -/
theorem sum_eq_zero_of_laplacian_mulVec_eq {a b : V → ℝ} (h : M.laplacian *ᵥ a = b) :
    ∑ v, b v = 0 :=
  Solver.sum_eq_zero_of_mulVec_eq M.laplacian_vecMul_one h

/-- **The pinned solve.** On a connected valid mesh, if `Σ b = 0` then for every vertex `v₀`
there is exactly one `a` with `a(v₀) = 0` and `L a = b`. It is obtained from the system without
the row and column of `v₀` (positive definite); the dropped equation holds automatically.

This is how the code solves singular Laplace systems: `DiskBff::solve_neumann` pins the last
boundary vertex, and `LaplaceSolver` in `crates/ddg-optimal-cones` drops the first unknown of a
closed mesh.

[formula: bff-neumann-to-dirichlet] -/
theorem existsUnique_pinned (hM : M.Valid) (hconn : M.Connected) {b : V → ℝ}
    (hb : ∑ v, b v = 0) (v₀ : V) : ∃! a : V → ℝ, a v₀ = 0 ∧ M.laplacian *ᵥ a = b := by
  obtain ⟨a, ⟨ha1, ha2⟩, hauniq⟩ :=
    M.existsUnique_poisson hM hconn (fun v => v ≠ v₀) ⟨v₀, by simp⟩ b (fun _ => 0)
  have h0 : a v₀ = 0 := ha1 v₀ (by simp)
  refine ⟨a, ⟨h0, ?_⟩, fun a' ha' => ?_⟩
  · funext v
    by_cases hv : v = v₀
    · subst hv
      exact Solver.mulVec_eq_of_sum_eq M.laplacian_vecMul_one hb v ha2
    · exact ha2 v hv
  · refine hauniq a' ⟨fun v hv => ?_, fun v _ => congrFun ha'.2 v⟩
    have : v = v₀ := by simpa using hv
    rw [this]; exact ha'.1

/-- **Existence for the Neumann problem**: on a connected valid mesh `L a = b` has a solution
iff `Σ b = 0`. -/
theorem exists_laplacian_mulVec_eq_iff (hM : M.Valid) (hconn : M.Connected) (b : V → ℝ) :
    (∃ a, M.laplacian *ᵥ a = b) ↔ ∑ v, b v = 0 := by
  constructor
  · rintro ⟨a, ha⟩; exact M.sum_eq_zero_of_laplacian_mulVec_eq ha
  · intro hb
    rcases isEmpty_or_nonempty V with hV | ⟨⟨v₀⟩⟩
    · exact ⟨fun _ => 0, funext fun v => (hV.false v).elim⟩
    · obtain ⟨a, ⟨_, ha⟩, _⟩ := M.existsUnique_pinned hM hconn hb v₀
      exact ⟨a, ha⟩

/-- Two solutions of `L a = b` on a connected valid mesh differ by a constant. -/
theorem sub_const_of_laplacian_mulVec_eq (hM : M.Valid) (hconn : M.Connected) {a a' b : V → ℝ}
    (h : M.laplacian *ᵥ a = b) (h' : M.laplacian *ᵥ a' = b) (i j : V) :
    a' i - a i = a' j - a j := by
  have : M.laplacian *ᵥ (a' - a) = 0 := by rw [Matrix.mulVec_sub, h, h']; simp
  exact (M.laplacian_mulVec_eq_zero_iff hM hconn _).mp this i j

/-- **The Neumann-to-Dirichlet map is well defined.** On a connected valid mesh, if `Σ b = 0`
and `m` are weights with `Σ m ≠ 0` (the lumped masses in `neumann_to_dirichlet`; `wE·m` for the
scale factor of a closed mesh), there is exactly one `a` with `L a = b` and `Σ_v m_v a_v = 0`.
For the Neumann problem of the code, `b = φ − [0; h]`, and `Σ b = 0` reads `Σ φ = Σ h`.

[formula: bff-neumann-to-dirichlet] -/
theorem existsUnique_neumann (hM : M.Valid) (hconn : M.Connected) {b : V → ℝ}
    (hb : ∑ v, b v = 0) (m : V → ℝ) (hm : ∑ v, m v ≠ 0) :
    ∃! a : V → ℝ, M.laplacian *ᵥ a = b ∧ ∑ v, m v * a v = 0 := by
  obtain ⟨a₀, ha₀⟩ := (M.exists_laplacian_mulVec_eq_iff hM hconn b).mpr hb
  set c : ℝ := -(∑ v, m v * a₀ v) / ∑ v, m v with hc
  refine ⟨fun v => a₀ v + c, ⟨?_, ?_⟩, ?_⟩
  · have : (fun v => a₀ v + c) = a₀ + c • (fun _ => (1 : ℝ)) := by ext v; simp
    rw [this, Matrix.mulVec_add, Matrix.mulVec_smul, M.laplacian_mulVec_one, ha₀]
    simp
  · simp only [mul_add, Finset.sum_add_distrib, ← Finset.sum_mul]
    rw [hc]; field_simp; ring
  · rintro a ⟨ha1, ha2⟩
    have hconst := M.sub_const_of_laplacian_mulVec_eq hM hconn ha₀ ha1
    funext v
    have h1 : ∀ j, a j = a₀ j + (a v - a₀ v) := fun j => by have := hconst v j; linarith
    have h2 : ∑ j, m j * a j = ∑ j, m j * a₀ j + (a v - a₀ v) * ∑ j, m j := by
      rw [Finset.mul_sum, ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun j _ => ?_
      rw [h1 j]; ring
    rw [ha2] at h2
    have h3 : a v - a₀ v = c := by
      rw [hc]; field_simp; linarith
    show a v = a₀ v + c
    linarith

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
