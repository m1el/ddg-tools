import DdgProofs.GaussBonnet

/-!
# Consistent orientation of the faces of a `TriMesh`

A `TriMesh` already lists the corners of every face in an order, `a, b, c`. Read as a cycle this
gives the face three **halfedges** `a → b`, `b → c`, `c → a`; the halfedge `(f, k)` runs from
corner `k` to corner `k + 1` (mod 3), the numbering `3f + k` of the Rust complexes
(`crates/ddg-intrinsic/src/complex.rs`, `crates/ddg-core`).

* `TriMesh.Oriented`: no directed edge is traversed by two halfedges (the map halfedge ↦
  (tail, head) is injective). Together with `TriMesh.Manifold` this is the usual notion: an edge
  between two faces is traversed once in each direction, a boundary edge once.
* `HasTwin h`: the reversed halfedge exists; `twin` returns it (and `h` itself on the boundary);
  on an oriented mesh `twin` is an involution (`twinPerm`).
* `not_hasTwin_iff`: on a valid oriented mesh a halfedge has no twin iff its edge is a boundary
  edge in the sense of `GaussBonnet.lean` (it lies in one face). Hence on a closed mesh every
  halfedge has a twin (`hasTwin_of_closed`).
* `sum_halfedge_antisymm`: for an antisymmetric `w`, `Σ_h w(tail h, head h)` is the sum over the
  boundary halfedges alone (the interior terms cancel in pairs). This is the telescoping behind
  the shoelace formula of `Bff/ConformalEnergy.lean` and the Poincaré–Hopf sum of
  `Intrinsic/SingularityIndex.lean`.

The existing theorems (Gauss–Bonnet, the cotan Laplacian) do not use the orientation; nothing in
`TriMesh` changes, `Oriented` is an extra hypothesis where it is needed.
-/

noncomputable section

namespace DdgProofs

namespace TriMesh

variable {V F : Type*} (M : TriMesh V F)

/-- The corners `a, b, c` of the face `f`, as a function of `Fin 3`. -/
def corner (f : F) : Fin 3 → V := ![M.a f, M.b f, M.c f]

/-- The halfedge `(f, k)` as a directed edge: from corner `k` to corner `k + 1` of `f`. -/
def halfedge (h : F × Fin 3) : V × V := (M.corner h.1 h.2, M.corner h.1 (h.2 + 1))

@[simp] theorem halfedge_zero (f : F) : M.halfedge (f, 0) = (M.a f, M.b f) := rfl
@[simp] theorem halfedge_one (f : F) : M.halfedge (f, 1) = (M.b f, M.c f) := rfl
@[simp] theorem halfedge_two (f : F) : M.halfedge (f, 2) = (M.c f, M.a f) := rfl

/-- **Consistent orientation**: no directed edge is traversed by two halfedges. -/
def Oriented : Prop := Function.Injective M.halfedge

/-- The halfedge `h` has a twin: some halfedge traverses its edge in the opposite direction. -/
def HasTwin (h : F × Fin 3) : Prop := ∃ h', M.halfedge h' = (M.halfedge h).swap

open Classical in
/-- The twin of `h`, and `h` itself if there is none (a boundary halfedge). -/
def twin (h : F × Fin 3) : F × Fin 3 := if ht : M.HasTwin h then ht.choose else h

theorem halfedge_twin {h : F × Fin 3} (ht : M.HasTwin h) :
    M.halfedge (M.twin h) = (M.halfedge h).swap := by
  rw [twin, dite_eq_left_of_eq_true (eq_true ht)]
  exact ht.choose_spec

theorem twin_of_not_hasTwin {h : F × Fin 3} (ht : ¬M.HasTwin h) : M.twin h = h := by
  rw [twin, dite_eq_right_of_eq_false (eq_false ht)]

theorem hasTwin_twin {h : F × Fin 3} (ht : M.HasTwin h) : M.HasTwin (M.twin h) :=
  ⟨h, by rw [M.halfedge_twin ht, Prod.swap_swap]⟩

theorem twin_twin (hO : M.Oriented) (h : F × Fin 3) : M.twin (M.twin h) = h := by
  by_cases ht : M.HasTwin h
  · apply hO
    rw [M.halfedge_twin (M.hasTwin_twin ht), M.halfedge_twin ht, Prod.swap_swap]
  · rw [M.twin_of_not_hasTwin ht, M.twin_of_not_hasTwin ht]

/-- On an oriented mesh `twin` is an involution of the halfedges. -/
def twinPerm (hO : M.Oriented) : Equiv.Perm (F × Fin 3) :=
  Function.Involutive.toPerm M.twin (M.twin_twin hO)

@[simp] theorem twinPerm_apply (hO : M.Oriented) (h : F × Fin 3) : M.twinPerm hO h = M.twin h :=
  rfl

/-! ## Halfedges and the unordered edges of `GaussBonnet.lean` -/

variable [DecidableEq V]

theorem mem_faceEdges_iff (f : F) (x y : V) :
    s(x, y) ∈ M.faceEdges f ↔ ∃ k, M.halfedge (f, k) = (x, y) ∨ M.halfedge (f, k) = (y, x) := by
  simp only [faceEdges, Finset.mem_insert, Finset.mem_singleton, Sym2.eq_iff]
  constructor
  · rintro ((⟨rfl, rfl⟩ | ⟨rfl, rfl⟩) | (⟨rfl, rfl⟩ | ⟨rfl, rfl⟩) | (⟨rfl, rfl⟩ | ⟨rfl, rfl⟩))
    · exact ⟨0, Or.inl rfl⟩
    · exact ⟨0, Or.inr rfl⟩
    · exact ⟨1, Or.inl rfl⟩
    · exact ⟨1, Or.inr rfl⟩
    · exact ⟨2, Or.inl rfl⟩
    · exact ⟨2, Or.inr rfl⟩
  · rintro ⟨k, hk⟩
    fin_cases k <;>
    simp only [Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk, halfedge_zero, halfedge_one,
      halfedge_two, Prod.mk.injEq] at hk <;>
    rcases hk with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;> simp

theorem halfedge_mem_faceEdges (h : F × Fin 3) :
    s((M.halfedge h).1, (M.halfedge h).2) ∈ M.faceEdges h.1 :=
  (M.mem_faceEdges_iff _ _ _).2 ⟨h.2, Or.inl rfl⟩

omit [DecidableEq V] in
/-- In a valid face the reversal of a halfedge is not a halfedge of the same face. -/
theorem halfedge_ne_swap (hM : M.Valid) (f : F) (k l : Fin 3) :
    M.halfedge (f, l) ≠ (M.halfedge (f, k)).swap := by
  have hab := hM.ab f
  have hbc := hM.bc f
  have hca := hM.ca f
  fin_cases k <;> fin_cases l <;>
  simp only [Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk, halfedge_zero, halfedge_one,
    halfedge_two, Prod.swap_prod_mk, ne_eq, Prod.mk.injEq, not_and] <;>
  intros <;> simp_all [eq_comm]

variable [Fintype F]

/-- **A halfedge has no twin iff its edge is a boundary edge** (lies in exactly one face), on a
valid oriented mesh. -/
theorem not_hasTwin_iff (hM : M.Valid) (hO : M.Oriented) (h : F × Fin 3) :
    ¬M.HasTwin h ↔ s((M.halfedge h).1, (M.halfedge h).2) ∈ M.boundaryEdges := by
  set x := (M.halfedge h).1
  set y := (M.halfedge h).2
  have hxy : M.halfedge h = (x, y) := rfl
  have hedge : s(x, y) ∈ M.edges := by
    simp only [edges, Finset.mem_biUnion, Finset.mem_univ, true_and]
    exact ⟨h.1, M.halfedge_mem_faceEdges h⟩
  simp only [boundaryEdges, Finset.mem_filter, hedge, true_and]
  unfold edgeDegree
  constructor
  · intro ht
    rw [Finset.card_eq_one]
    refine ⟨h.1, ?_⟩
    ext g
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_singleton]
    constructor
    · intro hg
      obtain ⟨l, hl | hl⟩ := (M.mem_faceEdges_iff g x y).1 hg
      · exact congrArg Prod.fst (hO (hl.trans hxy.symm))
      · exact (ht ⟨(g, l), by rw [hl, hxy]; rfl⟩).elim
    · rintro rfl
      exact M.halfedge_mem_faceEdges h
  · rintro hcard ⟨⟨g, l⟩, hgl⟩
    have hg : s(x, y) ∈ M.faceEdges g :=
      (M.mem_faceEdges_iff g x y).2 ⟨l, Or.inr (by rw [hgl, hxy]; rfl)⟩
    have hf : s(x, y) ∈ M.faceEdges h.1 := M.halfedge_mem_faceEdges h
    obtain ⟨f0, hf0⟩ := Finset.card_eq_one.1 hcard
    have e1 : g = f0 := by
      have : g ∈ Finset.univ.filter fun f => s(x, y) ∈ M.faceEdges f := by simp [hg]
      rw [hf0] at this
      simpa using this
    have e2 : h.1 = f0 := by
      have : h.1 ∈ Finset.univ.filter fun f => s(x, y) ∈ M.faceEdges f := by simp [hf]
      rw [hf0] at this
      simpa using this
    have : M.halfedge (h.1, l) = (M.halfedge h).swap := by
      rw [e2, ← e1]
      exact hgl
    exact M.halfedge_ne_swap hM h.1 h.2 l this

/-- **On a closed oriented mesh every halfedge has a twin.** -/
theorem hasTwin_of_closed [Fintype V] (hM : M.Valid) (hO : M.Oriented) (hclosed : M.boundaryVertices = ∅)
    (h : F × Fin 3) : M.HasTwin h := by
  by_contra ht
  have hb := (M.not_hasTwin_iff hM hO h).1 ht
  have : (M.halfedge h).1 ∈ M.boundaryVertices := by
    simp only [boundaryVertices, Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨_, hb, Sym2.mem_mk_left _ _⟩
  simp [hclosed] at this

omit [DecidableEq V] in
open Classical in
/-- **Telescoping over an oriented mesh**: for an antisymmetric `w` (`w y x = −w x y`) the terms
of two twin halfedges cancel, so the sum over all halfedges is the sum over the halfedges without
a twin (the boundary halfedges). -/
theorem sum_halfedge_antisymm (hO : M.Oriented) (w : V → V → ℝ)
    (hw : ∀ x y, w y x = -w x y) :
    ∑ h, w (M.halfedge h).1 (M.halfedge h).2 =
      ∑ h ∈ Finset.univ.filter (fun h => ¬M.HasTwin h), w (M.halfedge h).1 (M.halfedge h).2 := by
  classical
  set g : F × Fin 3 → ℝ := fun h => w (M.halfedge h).1 (M.halfedge h).2 with hg
  have hcomp : ∑ h, g (M.twinPerm hO h) = ∑ h, g h := Equiv.sum_comp _ g
  have htw : ∀ h, g (M.twinPerm hO h) = if M.HasTwin h then -g h else g h := by
    intro h
    split_ifs with ht
    · simp only [hg, twinPerm_apply, M.halfedge_twin ht, Prod.fst_swap, Prod.snd_swap]
      exact hw _ _
    · simp only [twinPerm_apply, M.twin_of_not_hasTwin ht]
  simp only [htw] at hcomp
  rw [Finset.sum_ite] at hcomp
  have hsplit := Finset.sum_filter_add_sum_filter_not Finset.univ (fun h => M.HasTwin h) g
  simp only [Finset.sum_neg_distrib] at hcomp
  show ∑ h, g h = ∑ h ∈ Finset.univ.filter (fun h => ¬M.HasTwin h), g h
  linarith

/-! ## Examples: the hypotheses are satisfiable -/

end TriMesh

/-- A consistently oriented regular tetrahedron: faces `012, 023, 031, 132`. -/
def exampleTetrahedronOriented : TriMesh (Fin 4) (Fin 4) :=
  ⟨![0, 0, 0, 1], ![1, 2, 3, 3], ![2, 3, 1, 2], fun _ _ => 1⟩

/-- It is oriented and closed, so every halfedge has a twin. -/
example : exampleTetrahedronOriented.Oriented ∧ exampleTetrahedronOriented.boundaryVertices = ∅ := by
  refine ⟨?_, by decide⟩
  intro h h' e
  revert h h'
  decide

/-- The tetrahedron of `GaussBonnet.lean` (faces `012, 023, 031, 123`) satisfies `Valid` and
`Manifold` but is **not** oriented: `012` and `123` both traverse `1 → 2`. Orientation is an extra
hypothesis. -/
example : ¬exampleTetrahedron.Oriented := by
  intro hO
  have := @hO (0, 1) (3, 0) (by decide)
  simp at this

/-- The single triangle is oriented, and none of its halfedges has a twin. -/
example : exampleTriangle.Oriented ∧ ∀ h, ¬exampleTriangle.HasTwin h := by
  refine ⟨?_, ?_⟩
  · intro h h' e
    revert h h'
    decide
  · intro h ⟨h', e⟩
    revert h h'
    decide

end DdgProofs
