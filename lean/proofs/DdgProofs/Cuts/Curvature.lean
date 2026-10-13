import DdgProofs.Oracle.Mesh

/-!
# Cut flow: integrated Gaussian curvature of a face

`crates/ddg-variational-cuts/SPEC.md` §3.1 (paper Sec. 5.2, Eq. 24 with the sign of D1). For a
triangle mesh with intrinsic edge lengths and a set `B` of vertices that carry no curvature (the
boundary vertices of `M`):

```
θ̃_i^f = θ_i^f · 2π/Θ_i   (i ∉ B),      θ̃_i^f = θ_i^f   (i ∈ B)
Ω_f   = θ̃_i^f + θ̃_j^f + θ̃_k^f − π
```

* **`cut-face-curvature`**: the rescaled angles around a vertex `i ∉ B` sum to `2π`
  (`TriMesh.sum_normAngleAt_of_not_mem`); around `i ∈ B` they sum to `Θ_i`; a face whose three
  vertices are in `B` has `Ω_f = 0`.
* **`cut-face-curvature-sum`**: `Σ_f Ω_f = Σ_{i ∉ B} (2π − Θ_i)` (`TriMesh.sum_faceCurvature`).
  With `B` the boundary vertices this is the sum of the angle defects of the interior vertices,
  `2πχ − Σ_B k` by Gauss–Bonnet, and `2πχ` on a closed mesh
  (`TriMesh.sum_faceCurvature_add_geodesicCurvature`, `TriMesh.sum_faceCurvature_closed`).
* `Cuts.face_omega_eq`, `Cuts.sum_face_omega`: the same about the oracle's `Ddg.Face.omega`.

**Hypothesis that SPEC does not state.** The identity needs `Θ_i ≠ 0` at every `i ∉ B`, i.e.
every such vertex lies in a face (`TriMesh.angleSum_pos`). A vertex in no face contributes `2π`
to the right-hand side and nothing to the left.
-/

noncomputable section

namespace DdgProofs

open Real Matrix

/-- A sum over all vertices of a function that vanishes off three distinct vertices. -/
theorem sum_eq_add_add_of_ne {V : Type*} [Fintype V] [DecidableEq V] {g : V → ℝ} {a b c : V}
    (hab : a ≠ b) (hbc : b ≠ c) (hca : c ≠ a)
    (h0 : ∀ v, v ≠ a → v ≠ b → v ≠ c → g v = 0) : ∑ v, g v = g a + g b + g c := by
  have h := Finset.sum_subset (Finset.subset_univ ({a, b, c} : Finset V)) (f := g)
    (fun v _ hv => by
      simp only [Finset.mem_insert, Finset.mem_singleton, not_or] at hv
      exact h0 v hv.1 hv.2.1 hv.2.2)
  rw [← h, Finset.sum_insert (by simp [hab, hca.symm]), Finset.sum_pair hbc]
  ring

namespace TriMesh

variable {V F : Type*} [DecidableEq V] [Fintype V] [Fintype F] (M : TriMesh V F)

omit [Fintype V] [Fintype F] in
theorem angleAt_nonneg (hM : M.Valid) (f : F) (v : V) : 0 ≤ M.angleAt f v := by
  unfold angleAt
  split_ifs
  · exact (cornerAngle_pos (hM.tri f)).le
  · exact (cornerAngle_pos (hM.tri f).rotate).le
  · exact (cornerAngle_pos (hM.tri f).rotate.rotate).le
  · exact le_rfl

omit [Fintype V] [Fintype F] in
theorem angleAt_pos (hM : M.Valid) {f : F} {v : V}
    (hv : v = M.a f ∨ v = M.b f ∨ v = M.c f) : 0 < M.angleAt f v := by
  unfold angleAt
  split_ifs with h1 h2 h3
  · exact cornerAngle_pos (hM.tri f)
  · exact cornerAngle_pos (hM.tri f).rotate
  · exact cornerAngle_pos (hM.tri f).rotate.rotate
  · exact absurd hv (by tauto)

omit [Fintype V] in
/-- The angle sum at a vertex of a face of a valid mesh is positive. -/
theorem angleSum_pos (hM : M.Valid) {v : V} (hv : ∃ f, v = M.a f ∨ v = M.b f ∨ v = M.c f) :
    0 < M.angleSum v := by
  obtain ⟨f, hf⟩ := hv
  exact Finset.sum_pos' (fun g _ => M.angleAt_nonneg hM g v)
    ⟨f, Finset.mem_univ f, M.angleAt_pos hM hf⟩

/-- Corner angle of `f` at `v`, rescaled so that the angles around a vertex outside `B` sum to
`2π`: `θ̃ = θ · 2π/Θ_v` for `v ∉ B`, `θ̃ = θ` for `v ∈ B` (SPEC §3.1). -/
def normAngleAt (B : Finset V) (f : F) (v : V) : ℝ :=
  if v ∈ B then M.angleAt f v else M.angleAt f v * (2 * π) / M.angleSum v

/-- Integrated Gaussian curvature of the face `f`: `Ω_f = θ̃_a + θ̃_b + θ̃_c − π` (SPEC §3.1). -/
def faceCurvature (B : Finset V) (f : F) : ℝ :=
  M.normAngleAt B f (M.a f) + M.normAngleAt B f (M.b f) + M.normAngleAt B f (M.c f) - π

omit [Fintype V] in
/-- **The rescaled angles around a vertex outside `B` sum to `2π`.**

[formula: cut-face-curvature] -/
theorem sum_normAngleAt_of_not_mem {B : Finset V} {v : V} (hv : v ∉ B)
    (hΘ : M.angleSum v ≠ 0) : ∑ f, M.normAngleAt B f v = 2 * π := by
  simp only [normAngleAt, ite_eq_right hv]
  rw [← Finset.sum_div, ← Finset.sum_mul]
  change M.angleSum v * (2 * π) / M.angleSum v = 2 * π
  field_simp

omit [Fintype V] in
/-- Around a vertex of `B` the angles are not rescaled: they sum to `Θ_v`. -/
theorem sum_normAngleAt_of_mem {B : Finset V} {v : V} (hv : v ∈ B) :
    ∑ f, M.normAngleAt B f v = M.angleSum v := by
  simp only [normAngleAt, ite_eq_left hv]
  rfl

omit [Fintype V] in
/-- The rescaled angle vanishes at a vertex that does not belong to the face. -/
theorem normAngleAt_of_ne (B : Finset V) {f : F} {v : V} (h1 : v ≠ M.a f) (h2 : v ≠ M.b f)
    (h3 : v ≠ M.c f) : M.normAngleAt B f v = 0 := by
  have : M.angleAt f v = 0 := by simp [angleAt, h1, h2, h3]
  simp [normAngleAt, this]

omit [Fintype V] in
/-- **A face whose three vertices are in `B` carries no curvature** (boundary vertices carry no
Gaussian curvature, SPEC D9a): its three angles are not rescaled and sum to `π`.

[formula: cut-face-curvature] -/
theorem faceCurvature_of_mem (hM : M.Valid) {B : Finset V} {f : F} (ha : M.a f ∈ B)
    (hb : M.b f ∈ B) (hc : M.c f ∈ B) : M.faceCurvature B f = 0 := by
  have h := cornerAngle_add_cornerAngle_add_cornerAngle (hM.tri f)
  have e1 : M.angleAt f (M.a f) = M.angleA f := by simp [angleAt]
  have e2 : M.angleAt f (M.b f) = M.angleB f := by simp [angleAt, (hM.ab f).symm]
  have e3 : M.angleAt f (M.c f) = M.angleC f := by simp [angleAt, hM.ca f, (hM.bc f).symm]
  simp only [faceCurvature, normAngleAt, ite_eq_left ha, ite_eq_left hb, ite_eq_left hc, e1, e2, e3]
  unfold angleA angleB angleC
  linarith

/-- **Total face curvature.** On a valid mesh, if every vertex outside `B` has a nonzero angle
sum (it lies in a face), then `Σ_f Ω_f = Σ_{v ∉ B} (2π − Θ_v)`.

[formula: cut-face-curvature-sum] -/
theorem sum_faceCurvature (hM : M.Valid) (B : Finset V)
    (hΘ : ∀ v, v ∉ B → M.angleSum v ≠ 0) :
    ∑ f, M.faceCurvature B f = ∑ v ∈ Bᶜ, (2 * π - M.angleSum v) := by
  have hface : ∀ f, M.faceCurvature B f = ∑ v, M.normAngleAt B f v - π := by
    intro f
    rw [sum_eq_add_add_of_ne (hM.ab f) (hM.bc f) (hM.ca f)
      (fun v h1 h2 h3 => M.normAngleAt_of_ne B h1 h2 h3)]
    rfl
  have hsum := M.sum_angleSum hM
  have h1 : ∑ f, M.faceCurvature B f =
      ∑ v, ∑ f, M.normAngleAt B f v - π * (Fintype.card F : ℝ) := by
    simp only [hface, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    rw [Finset.sum_comm]
    ring
  have h2 : ∑ v, ∑ f, M.normAngleAt B f v = ∑ v ∈ B, M.angleSum v + ∑ v ∈ Bᶜ, 2 * π := by
    rw [← Finset.sum_add_sum_compl B]
    congr 1
    · exact Finset.sum_congr rfl fun v hv => M.sum_normAngleAt_of_mem hv
    · exact Finset.sum_congr rfl fun v hv =>
        M.sum_normAngleAt_of_not_mem (Finset.mem_compl.mp hv) (hΘ v (Finset.mem_compl.mp hv))
  have h3 : ∑ v, M.angleSum v = ∑ v ∈ B, M.angleSum v + ∑ v ∈ Bᶜ, M.angleSum v :=
    (Finset.sum_add_sum_compl B _).symm
  rw [h1, h2, Finset.sum_sub_distrib]
  linarith

/-- With `B` the boundary vertices of the mesh, the total face curvature is the sum of the angle
defects of the interior vertices.

[formula: cut-face-curvature-sum] -/
theorem sum_faceCurvature_eq_sum_angleDefect (hM : M.Valid)
    (hΘ : ∀ v, v ∉ M.boundaryVertices → M.angleSum v ≠ 0) :
    ∑ f, M.faceCurvature M.boundaryVertices f = ∑ v ∈ M.boundaryVerticesᶜ, M.angleDefect v :=
  M.sum_faceCurvature hM _ hΘ

/-- **Gauss–Bonnet for the face curvatures**: on a valid manifold mesh whose interior vertices
all lie in a face, `Σ_f Ω_f + Σ_{v ∈ ∂M} k_v = 2πχ`.

[formula: cut-face-curvature-sum] -/
theorem sum_faceCurvature_add_geodesicCurvature (hM : M.Valid) (hman : M.Manifold)
    (hΘ : ∀ v, v ∉ M.boundaryVertices → M.angleSum v ≠ 0) :
    ∑ f, M.faceCurvature M.boundaryVertices f +
        ∑ v ∈ M.boundaryVertices, M.geodesicCurvature v = 2 * π * (M.eulerChar : ℝ) := by
  rw [M.sum_faceCurvature_eq_sum_angleDefect hM hΘ]
  exact M.gauss_bonnet hM hman

/-- **Closed mesh**: `Σ_f Ω_f = 2πχ`.

[formula: cut-face-curvature-sum] -/
theorem sum_faceCurvature_closed (hM : M.Valid) (hman : M.Manifold)
    (hclosed : M.boundaryVertices = ∅) (hΘ : ∀ v, M.angleSum v ≠ 0) :
    ∑ f, M.faceCurvature ∅ f = 2 * π * (M.eulerChar : ℝ) := by
  rw [M.sum_faceCurvature hM ∅ (fun v _ => hΘ v), ← M.gauss_bonnet_closed hM hman hclosed]
  simp [angleDefect]

end TriMesh

/-! ## The oracle's definitions -/

namespace Cuts

open scoped DdgProofs.Oracle

variable {faces : List Ddg.Face} {len : ℕ → ℕ → ℝ} {nv : ℕ}

/-- The oracle's `normalizedAngle` at `ℝ`: `θ · 2π / Θ`. -/
theorem normalizedAngle_eq (θ Θ : ℝ) : Ddg.normalizedAngle θ Θ = θ * (2 * π) / Θ := by
  simp only [Ddg.normalizedAngle, Oracle.two_eq, Oracle.pi_eq]

/-- The oracle's `faceCurvatureDensity`: `k_f · A_f = Ω_f` for a face of nonzero area. -/
theorem faceCurvatureDensity_mul_area (Ω : ℝ) {A : ℝ} (hA : A ≠ 0) :
    Ddg.faceCurvatureDensity Ω A * A = Ω := by
  unfold Ddg.faceCurvatureDensity
  field_simp

/-- **The oracle's rescaled angles around an interior vertex sum to `2π`**: for a vertex that
is not on the boundary, with `Θ` the sum of the angles `θ_k` around it.

[formula: cut-face-curvature] -/
theorem sum_cornerNormalized {ι : Type*} (s : Finset ι) (θ : ι → ℝ) (hΘ : ∑ k ∈ s, θ k ≠ 0) :
    ∑ k ∈ s, Ddg.cornerNormalized false (θ k) (∑ k ∈ s, θ k) = 2 * π := by
  have e : ∀ k, Ddg.cornerNormalized false (θ k) (∑ k ∈ s, θ k) =
      θ k * (2 * π) / ∑ k ∈ s, θ k := by
    intro k
    simp [Ddg.cornerNormalized, normalizedAngle_eq]
  simp only [e]
  rw [← Finset.sum_div, ← Finset.sum_mul]
  field_simp

theorem cornerNormalized_eq (hv : Oracle.MeshValid faces len nv) (B : Finset (Fin nv))
    (isB : ℕ → Bool) (hB : ∀ v : Fin nv, isB v = true ↔ v ∈ B) {i : ℕ} (hi : i < nv) (θ : ℝ) :
    Ddg.cornerNormalized (isB i) θ (Ddg.angleSum faces len i) =
      if (⟨i, hi⟩ : Fin nv) ∈ B then θ
      else θ * (2 * π) / (Oracle.toTriMesh faces len nv hv.lt).angleSum ⟨i, hi⟩ := by
  have h := Oracle.angleSum_eq hv ⟨i, hi⟩
  simp only at h
  have hB' := hB ⟨i, hi⟩
  simp only at hB'
  rw [h]
  unfold Ddg.cornerNormalized
  by_cases hm : (⟨i, hi⟩ : Fin nv) ∈ B
  · rw [ite_eq_left hm, ite_eq_left (hB'.mpr hm)]
  · have : ¬(isB i = true) := fun h' => hm (hB'.mp h')
    rw [ite_eq_right hm, ite_eq_right this, normalizedAngle_eq]

/-- **The oracle's `Face.omega` is `TriMesh.faceCurvature`** of the mesh it describes, when
`angleSum` is the oracle's angle sum and `isB` marks the set `B`. -/
theorem face_omega_eq (hv : Oracle.MeshValid faces len nv) (B : Finset (Fin nv))
    (isB : ℕ → Bool) (hB : ∀ v : Fin nv, isB v = true ↔ v ∈ B) (k : Fin faces.length) :
    (faces.get k).omega len (Ddg.angleSum faces len) isB =
      (Oracle.toTriMesh faces len nv hv.lt).faceCurvature B k := by
  have hM := Oracle.toTriMesh_valid hv
  have ht := hv.tri _ (List.get_mem faces k)
  have hlt := hv.lt _ (List.get_mem faces k)
  have ha : (Oracle.toTriMesh faces len nv hv.lt).angleAt k
      ((Oracle.toTriMesh faces len nv hv.lt).a k) =
      (Oracle.toTriMesh faces len nv hv.lt).angleA k := by simp [TriMesh.angleAt]
  have hb : (Oracle.toTriMesh faces len nv hv.lt).angleAt k
      ((Oracle.toTriMesh faces len nv hv.lt).b k) =
      (Oracle.toTriMesh faces len nv hv.lt).angleB k := by
    simp [TriMesh.angleAt, (hM.ab k).symm]
  have hc : (Oracle.toTriMesh faces len nv hv.lt).angleAt k
      ((Oracle.toTriMesh faces len nv hv.lt).c k) =
      (Oracle.toTriMesh faces len nv hv.lt).angleC k := by
    simp [TriMesh.angleAt, hM.ca k, (hM.bc k).symm]
  have e1 := cornerNormalized_eq hv B isB hB hlt.1 ((faces.get k).angleA len)
  have e2 := cornerNormalized_eq hv B isB hB hlt.2.1 ((faces.get k).angleB len)
  have e3 := cornerNormalized_eq hv B isB hB hlt.2.2 ((faces.get k).angleC len)
  unfold Ddg.Face.omega Ddg.faceCurvature
  rw [e1, e2, e3]
  unfold TriMesh.faceCurvature TriMesh.normAngleAt
  rw [ha, hb, hc]
  simp only [Oracle.pi_eq, Ddg.Face.angleA, Ddg.Face.angleB, Ddg.Face.angleC,
    Oracle.cornerAngle_eq ht, Oracle.cornerAngle_eq ht.rotate,
    Oracle.cornerAngle_eq ht.rotate.rotate]
  rfl

/-- **Total face curvature, on the oracle's definitions.** For a valid oracle mesh in which every
vertex outside `B` lies in a face,
`Σ_f Face.omega f = Σ_{v ∉ B} (2π − angleSum v)`.

[formula: cut-face-curvature-sum] -/
theorem sum_face_omega (hv : Oracle.MeshValid faces len nv) (B : Finset (Fin nv))
    (isB : ℕ → Bool) (hB : ∀ v : Fin nv, isB v = true ↔ v ∈ B)
    (hcov : ∀ v : Fin nv, v ∉ B → ∃ f ∈ faces, f.Has v) :
    Ddg.sumList (faces.map (Ddg.Face.omega len (Ddg.angleSum faces len) isB)) =
      ∑ v ∈ Bᶜ, (2 * π - Ddg.angleSum faces len v) := by
  rw [Oracle.sumList_map_eq]
  have hM := Oracle.toTriMesh_valid hv
  have hΘ : ∀ v, v ∉ B → (Oracle.toTriMesh faces len nv hv.lt).angleSum v ≠ 0 := by
    intro v hvB
    obtain ⟨f, hf, hfv⟩ := hcov v hvB
    obtain ⟨k, rfl⟩ := List.mem_iff_get.mp hf
    refine ((Oracle.toTriMesh faces len nv hv.lt).angleSum_pos hM ⟨k, ?_⟩).ne'
    rcases hfv with h | h | h
    · exact Or.inl (Fin.ext h)
    · exact Or.inr (Or.inl (Fin.ext h))
    · exact Or.inr (Or.inr (Fin.ext h))
  rw [Finset.sum_congr rfl fun k _ => face_omega_eq hv B isB hB k,
    (Oracle.toTriMesh faces len nv hv.lt).sum_faceCurvature hM B hΘ]
  exact Finset.sum_congr rfl fun v _ => by rw [Oracle.angleSum_eq hv v]

/-- **Closed mesh, on the oracle's definitions**: if the mesh is a closed combinatorial manifold
and every vertex lies in a face, `Σ_f Face.omega f = 2πχ` (no vertex is marked as boundary).

[formula: cut-face-curvature-sum] -/
theorem sum_face_omega_closed (hv : Oracle.MeshValid faces len nv)
    (hman : (Oracle.toTriMesh faces len nv hv.lt).Manifold)
    (hclosed : (Oracle.toTriMesh faces len nv hv.lt).boundaryVertices = ∅)
    (hcov : ∀ v : Fin nv, ∃ f ∈ faces, f.Has v) :
    Ddg.sumList (faces.map (Ddg.Face.omega len (Ddg.angleSum faces len) fun _ => false)) =
      2 * π * ((Oracle.toTriMesh faces len nv hv.lt).eulerChar : ℝ) := by
  rw [sum_face_omega hv ∅ (fun _ => false) (by simp) (fun v _ => hcov v),
    ← (Oracle.toTriMesh faces len nv hv.lt).gauss_bonnet_closed (Oracle.toTriMesh_valid hv) hman
      hclosed]
  simp only [Finset.compl_empty, TriMesh.angleDefect]
  exact Finset.sum_congr rfl fun v _ => by rw [Oracle.angleSum_eq hv v]

end Cuts

end DdgProofs
