import DDG.Mesh.Gluing

/-!
# The edge flip of gluing data

The flip of an edge of a triangulated surface (Sharp–Soliman–Crane 2019, Sec. 2.1 and Alg. 5;
Bobenko–Springborn 2007, Sec. 2; Indermitte–Liebling–Troyanov–Clémençon 2001, Sec. 4), on the
gluing data of `DDG/Mesh/Gluing.lean`, exactly as `Complex::flip_edge` of
`crates/ddg-intrinsic/src/complex.rs` does it.

**Setting.** `h₀ = (f₀, a)` runs `i → j` in `f₀ = (i, j, k)`, its twin `h₁ = (f₁, b)` runs `j → i`
in `f₁ = (j, i, l)`. Write `h₀ⁿ = next h₀` (`j → k`), `h₀ᵖ = prev h₀` (`k → i`), and likewise
`h₁ⁿ` (`i → l`), `h₁ᵖ` (`l → j`). After the flip `f₀ = (i, l, k)`, `f₁ = (j, k, l)`: the new edge
is `h₀ⁿ = l → k`, `h₁ⁿ = k → l`; the halfedge `i → l` now has the id `h₀` and `j → k` the id `h₁`;
`h₀ᵖ`, `h₁ᵖ` and all other halfedges keep their ids. In other words the new twin map is the
conjugate `σ ∘ twin ∘ σ⁻¹` of the old one by the 4-cycle `σ = (h₀ h₀ⁿ h₁ h₁ⁿ)` (`flipPerm`):
`σ` sends every surviving halfedge to its new id and the dying edge `{h₀, h₁}` to the new edge
`{h₀ⁿ, h₁ⁿ}` (`flip_twin_apply`). The Rust code computes the same map halfedge by halfedge
(`moved` is `σ` on the halfedges it is applied to; `flip_twin_self`, `flip_twin_twin`,
`flip_twin_next`, `flip_twin_twin_next`, `flip_twin_of_ne`), including the degenerate diamonds in
which the four outer sides are glued among themselves (e.g. `twin h₀ⁿ = h₁ⁿ`): conjugation needs no
case distinction.

**When the flip is defined** (`Flippable`): the two halfedges of the edge lie in distinct faces,
`(twin h₀).1 ≠ h₀.1`. This excludes boundary edges (`twin h₀ = h₀`) and an edge along which one face
is glued to itself, which happens exactly when an end point has a single corner (degree one;
Rust `FlipRefusal::BoundaryEdge` and `FlipRefusal::DegreeOneEndpoint`, paper Fig. 6 right). The
flip of such an edge is not a flip (there is no quadrilateral). Constrained edges are a Rust
policy, not topology. Everything else is allowed: `i = j`, `k = l`, `k ∈ {i, j}`, outer sides glued
to each other, self-loops and multiple edges before or after.

## Main results

* `flip_twin_apply`: `twin' (σ x) = σ (twin x)`; the flip is again gluing data (`flip`, an
  involution by construction) and flips the edge again (`flippable_flip`).
* `flipVertexEquiv`: the flip bijects the vertices. Every corner keeps its vertex: the old corner
  `x` goes to the new corner `flipCorner x` at the same vertex (`flipVertexEquiv_vertexOf`); every
  surviving halfedge keeps its tail and its head (`flip_vertexOf_surviving`); the four vertices of
  the diamond keep their classes and the new edge runs from `l` to `k` (`flip_diamond_vertices`).
  No hypothesis on the four vertices being distinct is needed.
* `flip_numVertices`, `flip_numEdges`, `Fintype.card F` unchanged, so `flip_eulerChar`: **χ is
  preserved**. `flip_numBoundary`, `flip_isBoundary_iff`, `flip_closed_iff`: **the boundary is
  unchanged** (as halfedges, relabelled by `σ`, with their tails and heads).
* Labellings (the Rust `faces` array): `flipLabel` is the Rust update (`flipLabel_eq`) and keeps a
  faithful onto labelling faithful and onto (`IsLabeling.flip`, `Faithful.flip`).
* **Orientation**: every gluing data is oriented by construction; a flip keeps every face
  positively oriented and every gluing orientation-reversing (`IsLabeling.twin` applied to
  `IsLabeling.flip`).
* `flipFlipIso`: **flipping the same edge twice** (the second time at its new halfedge `h₀ⁿ`, as
  the Rust edge id does) gives the original gluing data up to the explicit isomorphism that swaps
  `f₀` and `f₁` and rotates them; hence the same `V`, `E`, `F`, `χ`.
-/

namespace DDG.Delaunay

namespace Gluing

variable {F : Type*} [DecidableEq F] (G : Gluing F)

/-- The 4-cycle `h₀ ↦ next h₀ ↦ twin h₀ ↦ next (twin h₀) ↦ h₀`: old halfedge ids to new ones. -/
def flipPerm (h₀ : Halfedge F) : Equiv.Perm (Halfedge F) :=
  Equiv.swap h₀ (G.twin h₀).next * Equiv.swap h₀ (G.twin h₀) * Equiv.swap h₀ h₀.next

/-- **The flip** of the edge of `h₀`: the twin map conjugated by `flipPerm` (Rust
`Complex::flip_edge`). Only meaningful when `Flippable h₀`.

[formula: it-complex-flip] -/
def flip (h₀ : Halfedge F) : Gluing F where
  twin x := G.flipPerm h₀ (G.twin ((G.flipPerm h₀).symm x))
  twin_twin x := by simp

/-- The edge of `h₀` can be flipped: its two halfedges lie in distinct faces (so it is interior
and no face is glued to itself along it). -/
def Flippable (h₀ : Halfedge F) : Prop := (G.twin h₀).1 ≠ h₀.1

/-- The new twin map is the conjugate of the old one. -/
@[simp] theorem flip_twin_apply (h₀ x : Halfedge F) :
    (G.flip h₀).twin (G.flipPerm h₀ x) = G.flipPerm h₀ (G.twin x) := by
  simp [flip]

variable {G} {h₀ : Halfedge F}

namespace Flippable

omit [DecidableEq F] in
/-- A halfedge of `f₀` is not a halfedge of `f₁`. -/
theorem ne (hf : G.Flippable h₀) (x y : Halfedge F) (hx : x.1 = h₀.1 := by rfl)
    (hy : y.1 = (G.twin h₀).1 := by rfl) : x ≠ y := fun e => hf (by rw [← hy, ← e, hx])

omit [DecidableEq F] in
/-- A halfedge of `f₁` is not a halfedge of `f₀`. -/
theorem ne' (hf : G.Flippable h₀) (x y : Halfedge F) (hx : x.1 = (G.twin h₀).1 := by rfl)
    (hy : y.1 = h₀.1 := by rfl) : x ≠ y := (hf.ne y x hy hx).symm

omit [DecidableEq F] in
theorem twin_ne (hf : G.Flippable h₀) : G.twin h₀ ≠ h₀ := hf.ne' _ _

end Flippable

section Perm

variable (hf : G.Flippable h₀)
include hf

theorem flipPerm_self : G.flipPerm h₀ h₀ = h₀.next := by
  simp only [flipPerm, Equiv.Perm.mul_apply]
  rw [Equiv.swap_apply_left,
    Equiv.swap_apply_of_ne_of_ne (Halfedge.next_ne _) (hf.ne h₀.next (G.twin h₀)),
    Equiv.swap_apply_of_ne_of_ne (Halfedge.next_ne _) (hf.ne h₀.next (G.twin h₀).next)]

theorem flipPerm_next : G.flipPerm h₀ h₀.next = G.twin h₀ := by
  simp only [flipPerm, Equiv.Perm.mul_apply]
  rw [Equiv.swap_apply_right, Equiv.swap_apply_left,
    Equiv.swap_apply_of_ne_of_ne (hf.ne' (G.twin h₀) h₀) (Halfedge.ne_next _)]

theorem flipPerm_twin : G.flipPerm h₀ (G.twin h₀) = (G.twin h₀).next := by
  simp only [flipPerm, Equiv.Perm.mul_apply]
  rw [Equiv.swap_apply_of_ne_of_ne (hf.ne' (G.twin h₀) h₀) (hf.ne' (G.twin h₀) h₀.next),
    Equiv.swap_apply_right, Equiv.swap_apply_left]

theorem flipPerm_twin_next : G.flipPerm h₀ (G.twin h₀).next = h₀ := by
  simp only [flipPerm, Equiv.Perm.mul_apply]
  rw [Equiv.swap_apply_of_ne_of_ne (hf.ne' (G.twin h₀).next h₀)
      (hf.ne' (G.twin h₀).next h₀.next),
    Equiv.swap_apply_of_ne_of_ne (hf.ne' (G.twin h₀).next h₀) (Halfedge.next_ne _),
    Equiv.swap_apply_right]

omit hf in
theorem flipPerm_of_ne {x : Halfedge F} (x0 : x ≠ h₀) (x1 : x ≠ h₀.next) (x2 : x ≠ G.twin h₀)
    (x3 : x ≠ (G.twin h₀).next) : G.flipPerm h₀ x = x := by
  simp only [flipPerm, Equiv.Perm.mul_apply]
  rw [Equiv.swap_apply_of_ne_of_ne x0 x1, Equiv.swap_apply_of_ne_of_ne x0 x2,
    Equiv.swap_apply_of_ne_of_ne x0 x3]

theorem flipPerm_symm_next : (G.flipPerm h₀).symm h₀.next = h₀ := by
  rw [Equiv.symm_apply_eq, flipPerm_self hf]

theorem flipPerm_symm_twin : (G.flipPerm h₀).symm (G.twin h₀) = h₀.next := by
  rw [Equiv.symm_apply_eq, flipPerm_next hf]

/-! ### The new twins, halfedge by halfedge (as Rust computes them) -/

/-- The new edge: `l → k` and `k → l` are twins. -/
theorem flip_twin_next : (G.flip h₀).twin h₀.next = (G.twin h₀).next := by
  have := flip_twin_apply G h₀ h₀
  rwa [flipPerm_self hf, flipPerm_twin hf] at this

theorem flip_twin_twin_next : (G.flip h₀).twin (G.twin h₀).next = h₀.next := by
  have := flip_twin_apply G h₀ (G.twin h₀)
  rwa [flipPerm_twin hf, G.twin_twin, flipPerm_self hf] at this

/-- `i → l` (now `h₀`) is glued where `i → l` was (Rust `twin_il = moved(twin[h1n])`). -/
theorem flip_twin_self :
    (G.flip h₀).twin h₀ = G.flipPerm h₀ (G.twin (G.twin h₀).next) := by
  have := flip_twin_apply G h₀ (G.twin h₀).next
  rwa [flipPerm_twin_next hf] at this

/-- `j → k` (now `h₁`) is glued where `j → k` was (Rust `twin_jk = moved(twin[h0n])`). -/
theorem flip_twin_twin :
    (G.flip h₀).twin (G.twin h₀) = G.flipPerm h₀ (G.twin h₀.next) := by
  have := flip_twin_apply G h₀ h₀.next
  rwa [flipPerm_next hf] at this

omit hf in
/-- Every other halfedge keeps its twin, renamed by `σ` (Rust `twin[twin_il] = h0`, …). -/
theorem flip_twin_of_ne {x : Halfedge F} (x0 : x ≠ h₀) (x1 : x ≠ h₀.next) (x2 : x ≠ G.twin h₀)
    (x3 : x ≠ (G.twin h₀).next) : (G.flip h₀).twin x = G.flipPerm h₀ (G.twin x) := by
  conv_lhs => rw [← flipPerm_of_ne x0 x1 x2 x3]
  exact flip_twin_apply _ _ _

/-- The flipped edge can be flipped again, at its new halfedge `next h₀`. -/
theorem flippable_flip : (G.flip h₀).Flippable h₀.next := by
  rw [Flippable, flip_twin_next hf]
  exact hf

end Perm

/-! ## Vertices -/

/-- The corner of the flipped surface at the same vertex as the old corner `x`: the two corners
at `i` merge into `h₀`, the two at `j` into `twin h₀`; every other corner keeps its id. -/
def flipCorner (G : Gluing F) (h₀ x : Halfedge F) : Halfedge F :=
  if x = (G.twin h₀).next then h₀ else if x = h₀.next then G.twin h₀ else x

/-- The old corner at the same vertex as the new corner `y`: the new corners at `l` and `k`
(`next h₀`, `next (twin h₀)`) go to the old ones (`prev (twin h₀)`, `prev h₀`). -/
def flipCornerInv (G : Gluing F) (h₀ y : Halfedge F) : Halfedge F :=
  if y = h₀.next then (G.twin h₀).prev else if y = (G.twin h₀).next then h₀.prev else y

section Vertices

variable (hf : G.Flippable h₀)
include hf

omit hf in
theorem flipCorner_twin_next : G.flipCorner h₀ (G.twin h₀).next = h₀ := by simp [flipCorner]

theorem flipCorner_next : G.flipCorner h₀ h₀.next = G.twin h₀ := by
  simp [flipCorner, hf.ne h₀.next (G.twin h₀).next]

omit hf in
theorem flipCorner_of_ne {x : Halfedge F} (x1 : x ≠ (G.twin h₀).next) (x2 : x ≠ h₀.next) :
    G.flipCorner h₀ x = x := by
  simp [flipCorner, x1, x2]

omit hf in
theorem flipCornerInv_next : G.flipCornerInv h₀ h₀.next = (G.twin h₀).prev := by
  simp [flipCornerInv]

theorem flipCornerInv_twin_next : G.flipCornerInv h₀ (G.twin h₀).next = h₀.prev := by
  simp [flipCornerInv, hf.ne' (G.twin h₀).next h₀.next]

omit hf in
theorem flipCornerInv_of_ne {y : Halfedge F} (y1 : y ≠ h₀.next) (y2 : y ≠ (G.twin h₀).next) :
    G.flipCornerInv h₀ y = y := by
  simp [flipCornerInv, y1, y2]

omit hf in
/-- The generators of the new vertex relation, indexed by the old glued halfedges. -/
theorem flip_glued_iff (u y : Halfedge F) :
    (G.flip h₀).Glued (G.flipPerm h₀ u) y ↔
      G.twin u ≠ u ∧ y = (G.flipPerm h₀ (G.twin u)).next := by
  simp only [Glued, flip_twin_apply, ne_eq, (G.flipPerm h₀).injective.eq_iff]

omit hf in
theorem flip_glued_of (u : Halfedge F) (hu : G.twin u ≠ u) :
    (G.flip h₀).Glued (G.flipPerm h₀ u) (G.flipPerm h₀ (G.twin u)).next :=
  (flip_glued_iff u _).2 ⟨hu, rfl⟩

/-- New corners at `l`: `next h₀` (new `l → k`) follows `prev (twin h₀)`. -/
theorem flip_glued_l : (G.flip h₀).Glued h₀.next (G.twin h₀).prev := by
  have := flip_glued_of (G := G) (h₀ := h₀) h₀ hf.twin_ne
  rwa [flipPerm_self hf, flipPerm_twin hf, Halfedge.next_next] at this

/-- New corners at `k`: `next (twin h₀)` (new `k → l`) follows `prev h₀`. -/
theorem flip_glued_k : (G.flip h₀).Glued (G.twin h₀).next h₀.prev := by
  have hb : G.twin (G.twin h₀) ≠ G.twin h₀ := by rw [G.twin_twin]; exact hf.twin_ne.symm
  have := flip_glued_of (G := G) (h₀ := h₀) (G.twin h₀) hb
  rwa [flipPerm_twin hf, G.twin_twin, flipPerm_self hf, Halfedge.next_next] at this

omit [DecidableEq F] in
/-- Old corners at `i`. -/
theorem glued_i : G.Glued h₀ (G.twin h₀).next := ⟨hf.twin_ne, rfl⟩

omit [DecidableEq F] in
/-- Old corners at `j`. -/
theorem glued_j : G.Glued (G.twin h₀) h₀.next :=
  ⟨by rw [G.twin_twin]; exact hf.twin_ne.symm, by rw [G.twin_twin]⟩

open Relation in
/-- An old generator maps to a chain of new ones. -/
theorem flip_glued_corner {x y : Halfedge F} (h : G.Glued x y) :
    EqvGen (G.flip h₀).Glued (G.flipCorner h₀ x) (G.flipCorner h₀ y) := by
  obtain ⟨hx, rfl⟩ := h
  have N1 := EqvGen.rel _ _ (flip_glued_l hf)
  have N2 := EqvGen.rel _ _ (flip_glued_k hf)
  rcases eq_or_ne x h₀ with rfl | x0
  · rw [flipCorner_twin_next, flipCorner_of_ne (hf.ne x (G.twin x).next) (Halfedge.ne_next _)]
    exact EqvGen.refl _
  rcases eq_or_ne x (G.twin h₀) with rfl | x1
  · rw [G.twin_twin, flipCorner_next hf,
      flipCorner_of_ne (Halfedge.ne_next _) (hf.ne' (G.twin h₀) h₀.next)]
    exact EqvGen.refl _
  rcases eq_or_ne x h₀.next with rfl | x2
  · rw [flipCorner_next hf]
    have g := flip_glued_of (G := G) (h₀ := h₀) _ hx
    rw [flipPerm_next hf] at g
    rcases eq_or_ne (G.twin h₀.next) (G.twin h₀).next with hw | hw
    · rw [hw, Halfedge.next_next,
        flipCorner_of_ne (Halfedge.prev_ne_next _) (hf.ne' (G.twin h₀).prev h₀.next)]
      rw [hw, flipPerm_twin_next hf] at g
      exact EqvGen.trans _ _ _ (EqvGen.rel _ _ g) N1
    · have w0 : G.twin h₀.next ≠ h₀ := fun e =>
        absurd (G.twin_eq_iff.1 e) (hf.ne h₀.next (G.twin h₀))
      have w1 : G.twin h₀.next ≠ G.twin h₀ := fun e => Halfedge.next_ne h₀ (G.twin_injective e)
      rw [flipPerm_of_ne w0 hx w1 hw] at g
      rw [flipCorner_of_ne (fun e => w1 (Halfedge.next_injective e))
        (fun e => w0 (Halfedge.next_injective e))]
      exact EqvGen.rel _ _ g
  rcases eq_or_ne x (G.twin h₀).next with rfl | x3
  · rw [flipCorner_twin_next]
    have g := flip_glued_of (G := G) (h₀ := h₀) _ hx
    rw [flipPerm_twin_next hf] at g
    rcases eq_or_ne (G.twin (G.twin h₀).next) h₀.next with hw | hw
    · rw [hw, Halfedge.next_next,
        flipCorner_of_ne (hf.ne h₀.prev (G.twin h₀).next) (Halfedge.prev_ne_next _)]
      rw [hw, flipPerm_next hf] at g
      exact EqvGen.trans _ _ _ (EqvGen.rel _ _ g) N2
    · have w0 : G.twin (G.twin h₀).next ≠ h₀ := fun e =>
        Halfedge.next_ne (G.twin h₀) (G.twin_eq_iff.1 e)
      have w1 : G.twin (G.twin h₀).next ≠ G.twin h₀ := fun e =>
        absurd (G.twin_injective e) (hf.ne' (G.twin h₀).next h₀)
      rw [flipPerm_of_ne w0 hw w1 hx] at g
      rw [flipCorner_of_ne (fun e => w1 (Halfedge.next_injective e))
        (fun e => w0 (Halfedge.next_injective e))]
      exact EqvGen.rel _ _ g
  rw [flipCorner_of_ne x3 x2]
  have g := flip_glued_of (G := G) (h₀ := h₀) _ hx
  rw [flipPerm_of_ne x0 x2 x1 x3] at g
  rcases eq_or_ne (G.twin x) h₀.next with hw | hw
  · rw [hw, Halfedge.next_next,
      flipCorner_of_ne (hf.ne h₀.prev (G.twin h₀).next) (Halfedge.prev_ne_next _)]
    rw [hw, flipPerm_next hf] at g
    exact EqvGen.trans _ _ _ (EqvGen.rel _ _ g) N2
  rcases eq_or_ne (G.twin x) (G.twin h₀).next with hw' | hw'
  · rw [hw', Halfedge.next_next,
      flipCorner_of_ne (Halfedge.prev_ne_next _) (hf.ne' (G.twin h₀).prev h₀.next)]
    rw [hw', flipPerm_twin_next hf] at g
    exact EqvGen.trans _ _ _ (EqvGen.rel _ _ g) N1
  have w0 : G.twin x ≠ h₀ := fun e => x1 (G.twin_eq_iff.1 e)
  have w1 : G.twin x ≠ G.twin h₀ := fun e => x0 (G.twin_injective e)
  rw [flipPerm_of_ne w0 hw w1 hw'] at g
  rw [flipCorner_of_ne (fun e => w1 (Halfedge.next_injective e))
    (fun e => w0 (Halfedge.next_injective e))]
  exact EqvGen.rel _ _ g

open Relation in
/-- A new generator maps back to a chain of old ones. -/
theorem flip_glued_cornerInv {x y : Halfedge F} (h : (G.flip h₀).Glued x y) :
    EqvGen G.Glued (G.flipCornerInv h₀ x) (G.flipCornerInv h₀ y) := by
  obtain ⟨u, rfl⟩ := (G.flipPerm h₀).surjective x
  obtain ⟨hu, rfl⟩ := (flip_glued_iff u y).1 h
  have O1 := EqvGen.rel _ _ (glued_i hf)
  have O2 := EqvGen.rel _ _ (glued_j hf)
  rcases eq_or_ne u h₀ with rfl | u0
  · rw [flipPerm_self hf, flipCornerInv_next, flipPerm_twin hf, Halfedge.next_next,
      flipCornerInv_of_ne (hf.ne' (G.twin u).prev u.next) (Halfedge.prev_ne_next _)]
    exact EqvGen.refl _
  rcases eq_or_ne u (G.twin h₀) with rfl | u1
  · rw [flipPerm_twin hf, flipCornerInv_twin_next hf, G.twin_twin, flipPerm_self hf,
      Halfedge.next_next,
      flipCornerInv_of_ne (Halfedge.prev_ne_next _) (hf.ne h₀.prev (G.twin h₀).next)]
    exact EqvGen.refl _
  rcases eq_or_ne u h₀.next with rfl | u2
  · rw [flipPerm_next hf,
      flipCornerInv_of_ne (hf.ne' (G.twin h₀) h₀.next) (Halfedge.ne_next _)]
    rcases eq_or_ne (G.twin h₀.next) (G.twin h₀).next with hw | hw
    · rw [hw, flipPerm_twin_next hf, flipCornerInv_next]
      have g : G.Glued h₀.next (G.twin h₀).prev := ⟨hu, by rw [hw, Halfedge.next_next]⟩
      exact EqvGen.trans _ _ _ O2 (EqvGen.rel _ _ g)
    · have w0 : G.twin h₀.next ≠ h₀ := fun e =>
        absurd (G.twin_eq_iff.1 e) (hf.ne h₀.next (G.twin h₀))
      have w1 : G.twin h₀.next ≠ G.twin h₀ := fun e => Halfedge.next_ne h₀ (G.twin_injective e)
      rw [flipPerm_of_ne w0 hu w1 hw, flipCornerInv_of_ne
        (fun e => w0 (Halfedge.next_injective e)) (fun e => w1 (Halfedge.next_injective e))]
      exact EqvGen.trans _ _ _ O2 (EqvGen.rel _ _ ⟨hu, rfl⟩)
  rcases eq_or_ne u (G.twin h₀).next with rfl | u3
  · rw [flipPerm_twin_next hf,
      flipCornerInv_of_ne (Halfedge.ne_next _) (hf.ne h₀ (G.twin h₀).next)]
    rcases eq_or_ne (G.twin (G.twin h₀).next) h₀.next with hw | hw
    · rw [hw, flipPerm_next hf, flipCornerInv_twin_next hf]
      have g : G.Glued (G.twin h₀).next h₀.prev := ⟨hu, by rw [hw, Halfedge.next_next]⟩
      exact EqvGen.trans _ _ _ O1 (EqvGen.rel _ _ g)
    · have w0 : G.twin (G.twin h₀).next ≠ h₀ := fun e =>
        Halfedge.next_ne (G.twin h₀) (G.twin_eq_iff.1 e)
      have w1 : G.twin (G.twin h₀).next ≠ G.twin h₀ := fun e =>
        absurd (G.twin_injective e) (hf.ne' (G.twin h₀).next h₀)
      rw [flipPerm_of_ne w0 hw w1 hu, flipCornerInv_of_ne
        (fun e => w0 (Halfedge.next_injective e)) (fun e => w1 (Halfedge.next_injective e))]
      exact EqvGen.trans _ _ _ O1 (EqvGen.rel _ _ ⟨hu, rfl⟩)
  rw [flipPerm_of_ne u0 u2 u1 u3, flipCornerInv_of_ne u2 u3]
  rcases eq_or_ne (G.twin u) h₀.next with hw | hw
  · rw [hw, flipPerm_next hf, flipCornerInv_twin_next hf]
    exact EqvGen.rel _ _ ⟨hu, by rw [hw, Halfedge.next_next]⟩
  rcases eq_or_ne (G.twin u) (G.twin h₀).next with hw' | hw'
  · rw [hw', flipPerm_twin_next hf, flipCornerInv_next]
    exact EqvGen.rel _ _ ⟨hu, by rw [hw', Halfedge.next_next]⟩
  have w0 : G.twin u ≠ h₀ := fun e => u1 (G.twin_eq_iff.1 e)
  have w1 : G.twin u ≠ G.twin h₀ := fun e => u0 (G.twin_injective e)
  rw [flipPerm_of_ne w0 hw w1 hw', flipCornerInv_of_ne
    (fun e => w0 (Halfedge.next_injective e)) (fun e => w1 (Halfedge.next_injective e))]
  exact EqvGen.rel _ _ ⟨hu, rfl⟩

open Relation in
theorem flipCornerInv_flipCorner (x : Halfedge F) :
    EqvGen G.Glued (G.flipCornerInv h₀ (G.flipCorner h₀ x)) x := by
  rcases eq_or_ne x (G.twin h₀).next with rfl | x1
  · rw [flipCorner_twin_next, flipCornerInv_of_ne (Halfedge.ne_next _) (hf.ne h₀ (G.twin h₀).next)]
    exact EqvGen.rel _ _ (glued_i hf)
  rcases eq_or_ne x h₀.next with rfl | x2
  · rw [flipCorner_next hf, flipCornerInv_of_ne (hf.ne' (G.twin h₀) h₀.next) (Halfedge.ne_next _)]
    exact EqvGen.rel _ _ (glued_j hf)
  rw [flipCorner_of_ne x1 x2, flipCornerInv_of_ne x2 x1]
  exact EqvGen.refl _

open Relation in
theorem flipCorner_flipCornerInv (y : Halfedge F) :
    EqvGen (G.flip h₀).Glued (G.flipCorner h₀ (G.flipCornerInv h₀ y)) y := by
  rcases eq_or_ne y h₀.next with rfl | y1
  · rw [flipCornerInv_next,
      flipCorner_of_ne (Halfedge.prev_ne_next _) (hf.ne' (G.twin h₀).prev h₀.next)]
    exact EqvGen.symm _ _ (EqvGen.rel _ _ (flip_glued_l hf))
  rcases eq_or_ne y (G.twin h₀).next with rfl | y2
  · rw [flipCornerInv_twin_next hf,
      flipCorner_of_ne (hf.ne h₀.prev (G.twin h₀).next) (Halfedge.prev_ne_next _)]
    exact EqvGen.symm _ _ (EqvGen.rel _ _ (flip_glued_k hf))
  rw [flipCornerInv_of_ne y1 y2, flipCorner_of_ne y2 y1]
  exact EqvGen.refl _

/-- **The flip bijects the vertices**: a corner goes to the corner at the same vertex. -/
def flipVertexEquiv : G.Vertex ≃ (G.flip h₀).Vertex :=
  quotientEquivOfMaps _ _ (G.flipCorner h₀) (G.flipCornerInv h₀)
    (fun _ _ h => eqvGen_map (G.flipCorner h₀) (fun _ _ hr => flip_glued_corner hf hr) h)
    (fun _ _ h => eqvGen_map (G.flipCornerInv h₀) (fun _ _ hr => flip_glued_cornerInv hf hr) h)
    (flipCornerInv_flipCorner hf) (flipCorner_flipCornerInv hf)

theorem flipVertexEquiv_vertexOf (x : Halfedge F) :
    flipVertexEquiv hf (G.vertexOf x) = (G.flip h₀).vertexOf (G.flipCorner h₀ x) := rfl

theorem flipVertexEquiv_symm_vertexOf (y : Halfedge F) :
    (flipVertexEquiv hf).symm ((G.flip h₀).vertexOf y) = G.vertexOf (G.flipCornerInv h₀ y) :=
  rfl

/-- **Every surviving halfedge keeps its tail and its head** (renamed by `σ`): for `x` not on the
flipped edge, the vertex of `σ x` and of `next (σ x)` in the new surface are those of `x` and of
`next x`. This includes the four outer sides of the diamond and every boundary halfedge.

[formula: it-complex-flip] -/
theorem flip_vertexOf_surviving {x : Halfedge F} (x0 : x ≠ h₀) (x1 : x ≠ G.twin h₀) :
    flipVertexEquiv hf (G.vertexOf x) = (G.flip h₀).vertexOf (G.flipPerm h₀ x) ∧
      flipVertexEquiv hf (G.vertexOf x.next) =
        (G.flip h₀).vertexOf (G.flipPerm h₀ x).next := by
  rw [flipVertexEquiv_vertexOf, flipVertexEquiv_vertexOf]
  rcases eq_or_ne x h₀.next with rfl | x2
  · rw [flipCorner_next hf, flipPerm_next hf, Halfedge.next_next,
      flipCorner_of_ne (hf.ne h₀.prev (G.twin h₀).next) (Halfedge.prev_ne_next _)]
    exact ⟨rfl, ((G.flip h₀).vertexOf_eq_of_glued (flip_glued_k hf)).symm⟩
  rcases eq_or_ne x (G.twin h₀).next with rfl | x3
  · rw [flipCorner_twin_next, flipPerm_twin_next hf, Halfedge.next_next,
      flipCorner_of_ne (Halfedge.prev_ne_next _) (hf.ne' (G.twin h₀).prev h₀.next)]
    exact ⟨rfl, ((G.flip h₀).vertexOf_eq_of_glued (flip_glued_l hf)).symm⟩
  rw [flipCorner_of_ne x3 x2, flipPerm_of_ne x0 x2 x1 x3,
    flipCorner_of_ne (fun e => x1 (Halfedge.next_injective e))
      (fun e => x0 (Halfedge.next_injective e))]
  exact ⟨rfl, rfl⟩

/-- **The four vertices of the diamond keep their classes**, `i, j, k, l` (the tails of `h₀`,
`twin h₀`, `prev h₀`, `prev (twin h₀)`), whether or not they coincide, and **the new edge runs from
`l` to `k`**.

[formula: it-complex-flip] -/
theorem flip_diamond_vertices :
    flipVertexEquiv hf (G.vertexOf h₀) = (G.flip h₀).vertexOf h₀ ∧
      flipVertexEquiv hf (G.vertexOf (G.twin h₀)) = (G.flip h₀).vertexOf (G.twin h₀) ∧
      flipVertexEquiv hf (G.vertexOf h₀.prev) = (G.flip h₀).vertexOf h₀.prev ∧
      flipVertexEquiv hf (G.vertexOf (G.twin h₀).prev) =
        (G.flip h₀).vertexOf (G.twin h₀).prev ∧
      (G.flip h₀).vertexOf h₀.next = flipVertexEquiv hf (G.vertexOf (G.twin h₀).prev) ∧
      (G.flip h₀).vertexOf h₀.next.next = flipVertexEquiv hf (G.vertexOf h₀.prev) := by
  simp only [flipVertexEquiv_vertexOf]
  rw [flipCorner_of_ne (hf.ne h₀ (G.twin h₀).next) (Halfedge.ne_next _),
    flipCorner_of_ne (Halfedge.ne_next _) (hf.ne' (G.twin h₀) h₀.next),
    flipCorner_of_ne (hf.ne h₀.prev (G.twin h₀).next) (Halfedge.prev_ne_next _),
    flipCorner_of_ne (Halfedge.prev_ne_next _) (hf.ne' (G.twin h₀).prev h₀.next),
    Halfedge.next_next]
  exact ⟨rfl, rfl, rfl, rfl, (G.flip h₀).vertexOf_eq_of_glued (flip_glued_l hf), rfl⟩

/-- The flip keeps the number of vertices. -/
theorem flip_numVertices : (G.flip h₀).numVertices = G.numVertices :=
  (Nat.card_congr (flipVertexEquiv hf)).symm

end Vertices

/-! ## Edges, boundary, Euler characteristic -/

/-- The flip bijects the edges by `σ`; the flipped edge goes to the new edge. -/
def flipEdgeEquiv (G : Gluing F) (h₀ : Halfedge F) : G.Edge ≃ (G.flip h₀).Edge :=
  quotientEquivOfMaps _ _ (G.flipPerm h₀) (G.flipPerm h₀).symm
    (by
      rintro x y (rfl | rfl)
      · exact Or.inl rfl
      · exact Or.inr (flip_twin_apply _ _ _).symm)
    (by
      rintro x y (rfl | rfl)
      · exact Or.inl rfl
      · right
        obtain ⟨u, rfl⟩ := (G.flipPerm h₀).surjective x
        simp)
    (fun x => by simp only [Equiv.symm_apply_apply]; exact Or.inl rfl)
    (fun y => by simp only [Equiv.apply_symm_apply]; exact Or.inl rfl)

/-- The flip keeps the number of edges. -/
theorem flip_numEdges (h₀ : Halfedge F) : (G.flip h₀).numEdges = G.numEdges :=
  (Nat.card_congr (G.flipEdgeEquiv h₀)).symm

/-- **A flip keeps the Euler characteristic** `χ = V − E + F`.

[formula: it-complex-flip] [formula: it-complex-euler] -/
theorem flip_eulerChar [Fintype F] (hf : G.Flippable h₀) :
    (G.flip h₀).eulerChar = G.eulerChar := by
  simp only [eulerChar, flip_numVertices hf, flip_numEdges]

theorem flip_isBoundary_iff (h₀ x : Halfedge F) :
    (G.flip h₀).twin (G.flipPerm h₀ x) = G.flipPerm h₀ x ↔ G.twin x = x := by
  rw [flip_twin_apply, (G.flipPerm h₀).injective.eq_iff]

/-- **The boundary is unchanged**: `σ` maps the boundary halfedges onto the boundary halfedges
(each with its tail and head, `flip_vertexOf_surviving`). -/
theorem flip_numBoundary [Fintype F] (h₀ : Halfedge F) :
    (G.flip h₀).numBoundary = G.numBoundary := by
  symm
  apply Finset.card_nbij' (G.flipPerm h₀) (G.flipPerm h₀).symm
  · intro x hx
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at hx ⊢
    exact (flip_isBoundary_iff h₀ x).2 hx
  · intro y hy
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at hy ⊢
    rw [← flip_isBoundary_iff h₀, Equiv.apply_symm_apply]
    exact hy
  · intro x _; simp
  · intro y _; simp

theorem flip_closed_iff (h₀ : Halfedge F) : (G.flip h₀).Closed ↔ G.Closed := by
  constructor
  · intro hc x hx
    exact hc _ ((flip_isBoundary_iff h₀ x).2 hx)
  · intro hc y hy
    obtain ⟨x, rfl⟩ := (G.flipPerm h₀).surjective y
    exact hc x ((flip_isBoundary_iff h₀ x).1 hy)

/-! ## Labellings: the Rust `faces` array -/

/-- The labelling after the flip: `lab ∘ flipCornerInv`. -/
def flipLabel {V : Type*} (G : Gluing F) (h₀ : Halfedge F) (lab : Halfedge F → V) :
    Halfedge F → V :=
  fun y => lab (G.flipCornerInv h₀ y)

/-- `flipLabel` is the Rust update `faces[f0][(a+1)%3] = l; faces[f1][(b+1)%3] = k`: the corner
`next h₀` gets the label `l` of `prev (twin h₀)`, the corner `next (twin h₀)` the label `k` of
`prev h₀`, and every other corner keeps its label.

[formula: it-complex-flip] -/
theorem flipLabel_eq {V : Type*} (hf : G.Flippable h₀) (lab : Halfedge F → V) :
    G.flipLabel h₀ lab h₀.next = lab (G.twin h₀).prev ∧
      G.flipLabel h₀ lab (G.twin h₀).next = lab h₀.prev ∧
      ∀ y, y ≠ h₀.next → y ≠ (G.twin h₀).next → G.flipLabel h₀ lab y = lab y :=
  ⟨congrArg lab flipCornerInv_next, congrArg lab (flipCornerInv_twin_next hf),
    fun _ y1 y2 => congrArg lab (flipCornerInv_of_ne y1 y2)⟩

/-- A labelling stays a labelling. -/
theorem IsLabeling.flip {V : Type*} {lab : Halfedge F → V} (hl : G.IsLabeling lab)
    (hf : G.Flippable h₀) : (G.flip h₀).IsLabeling (G.flipLabel h₀ lab) :=
  fun _ _ h => hl.eq_of_eqvGen (flip_glued_cornerInv hf h)

/-- A faithful labelling stays faithful (each label is still one fan). -/
theorem Faithful.flip {V : Type*} {lab : Halfedge F → V} (hfa : G.Faithful lab)
    (hf : G.Flippable h₀) :
    (G.flip h₀).Faithful (G.flipLabel h₀ lab) := by
  intro x y hxy
  have h1 := G.vertexOf_eq_vertexOf.1 (hfa _ _ hxy)
  have h2 := eqvGen_map (G.flipCorner h₀) (fun _ _ hr => flip_glued_corner hf hr) h1
  have := Relation.EqvGen.trans _ _ _ (Relation.EqvGen.symm _ _ (flipCorner_flipCornerInv hf x))
    (Relation.EqvGen.trans _ _ _ h2 (flipCorner_flipCornerInv hf y))
  exact (G.flip h₀).vertexOf_eq_vertexOf.2 this

/-- An onto labelling stays onto. -/
theorem flipLabel_surjective {V : Type*} {lab : Halfedge F → V} (hl : G.IsLabeling lab)
    (hs : Function.Surjective lab) (hf : G.Flippable h₀) :
    Function.Surjective (G.flipLabel h₀ lab) := by
  intro v
  obtain ⟨x, rfl⟩ := hs v
  refine ⟨G.flipCorner h₀ x, ?_⟩
  exact hl.eq_of_eqvGen (flipCornerInv_flipCorner hf x)

/-! ## Flipping twice -/

/-- The isomorphism from `G` to the twice flipped `G`: faces `f₀ = h₀.1` and `f₁ = (twin h₀).1`
swapped, `f₀` rotated so that `h₀` lands on `prev (twin h₀)`, `f₁` so that `twin h₀` lands on
`prev h₀`. -/
def flipFlipShift (G : Gluing F) (h₀ : Halfedge F) (f : F) : Fin 3 :=
  if f = h₀.1 then (G.twin h₀).2 + 2 - h₀.2
  else if f = (G.twin h₀).1 then h₀.2 + 2 - (G.twin h₀).2 else 0

omit [DecidableEq F] in
private theorem fin3_aux (a b : Fin 3) :
    a + (b + 2 - a) = b + 2 ∧ a + 1 + (b + 2 - a) = b ∧ a + 2 + (b + 2 - a) = b + 1 := by
  revert a b; decide

/-- The relabelling of `flipFlipIso`, as a map of halfedges. -/
private def psi (G : Gluing F) (h₀ x : Halfedge F) : Halfedge F :=
  (Equiv.swap h₀.1 (G.twin h₀).1 x.1, x.2 + G.flipFlipShift h₀ x.1)

private theorem psi_other {x : Halfedge F} (h0 : x.1 ≠ h₀.1) (h1 : x.1 ≠ (G.twin h₀).1) :
    psi G h₀ x = x := by
  simp only [psi, flipFlipShift, h0, h1, ↓reduceIte, add_zero,
    Equiv.swap_apply_of_ne_of_ne h0 h1]

private theorem psi_f0 (x : Halfedge F) (hx : x.1 = h₀.1) :
    psi G h₀ x = ((G.twin h₀).1, x.2 + ((G.twin h₀).2 + 2 - h₀.2)) := by
  simp only [psi, hx, Equiv.swap_apply_left, flipFlipShift, ↓reduceIte]

private theorem psi_f1 (hf : G.Flippable h₀) (x : Halfedge F) (hx : x.1 = (G.twin h₀).1) :
    psi G h₀ x = (h₀.1, x.2 + (h₀.2 + 2 - (G.twin h₀).2)) := by
  simp only [psi, hx, Equiv.swap_apply_right, flipFlipShift,
    show (G.twin h₀).1 ≠ h₀.1 from hf, ↓reduceIte]

private theorem swap_twin_swap (x : Halfedge F) :
    Equiv.swap h₀ (G.twin h₀) (G.twin (Equiv.swap h₀ (G.twin h₀) x)) = G.twin x := by
  rcases eq_or_ne x h₀ with rfl | x0
  · rw [Equiv.swap_apply_left, G.twin_twin, Equiv.swap_apply_left]
  rcases eq_or_ne x (G.twin h₀) with rfl | x1
  · rw [Equiv.swap_apply_right, Equiv.swap_apply_right, G.twin_twin]
  rw [Equiv.swap_apply_of_ne_of_ne x0 x1]
  apply Equiv.swap_apply_of_ne_of_ne
  · exact fun e => x1 (G.twin_eq_iff.1 e)
  · exact fun e => x0 (G.twin_injective e)

section Twice

variable (hf : G.Flippable h₀)
include hf

/-- The composite of the two relabellings `σ₂ ∘ σ₁` is `ψ ∘ swap h₀ (twin h₀)`. -/
private theorem flipPerm_flipPerm (x : Halfedge F) :
    (G.flip h₀).flipPerm h₀.next (G.flipPerm h₀ x) =
      psi G h₀ (Equiv.swap h₀ (G.twin h₀) x) := by
  have hf2 := flippable_flip hf
  have t2 : (G.flip h₀).twin h₀.next = (G.twin h₀).next := flip_twin_next hf
  -- the second 4-cycle: next h₀ ↦ prev h₀ ↦ next (twin h₀) ↦ prev (twin h₀) ↦ next h₀
  have s0 : (G.flip h₀).flipPerm h₀.next h₀.next = h₀.prev := by
    rw [flipPerm_self hf2, Halfedge.next_next]
  have s1 : (G.flip h₀).flipPerm h₀.next h₀.prev = (G.twin h₀).next := by
    rw [← Halfedge.next_next, flipPerm_next hf2, t2]
  have s2 : (G.flip h₀).flipPerm h₀.next (G.twin h₀).next = (G.twin h₀).prev := by
    rw [← t2, flipPerm_twin hf2, t2, Halfedge.next_next]
  have s3 : (G.flip h₀).flipPerm h₀.next (G.twin h₀).prev = h₀.next := by
    rw [← Halfedge.next_next, ← t2, flipPerm_twin_next hf2]
  have s_fix : ∀ y : Halfedge F, y ≠ h₀.next → y ≠ h₀.prev → y ≠ (G.twin h₀).next →
      y ≠ (G.twin h₀).prev → (G.flip h₀).flipPerm h₀.next y = y := by
    intro y y0 y1 y2 y3
    apply flipPerm_of_ne y0 (by rwa [Halfedge.next_next]) (by rwa [t2])
      (by rwa [t2, Halfedge.next_next])
  have fx := fin3_aux h₀.2 (G.twin h₀).2
  have fy := fin3_aux (G.twin h₀).2 h₀.2
  by_cases x0 : x.1 = h₀.1
  · rcases Halfedge.eq_or_eq_next_or_eq_prev x0 with rfl | rfl | rfl
    · rw [Equiv.swap_apply_left, flipPerm_self hf, s0, psi_f1 hf (G.twin x) rfl]
      exact Prod.ext rfl fy.1.symm
    · rw [Equiv.swap_apply_of_ne_of_ne (Halfedge.next_ne _) (hf.ne h₀.next (G.twin h₀)),
        flipPerm_next hf, psi_f0 (h₀ := h₀) h₀.next rfl,
        s_fix _ (hf.ne' (G.twin h₀) h₀.next) (hf.ne' (G.twin h₀) h₀.prev) (Halfedge.ne_next _)
          (Halfedge.ne_prev _)]
      exact Prod.ext rfl fx.2.1.symm
    · rw [Equiv.swap_apply_of_ne_of_ne (Halfedge.prev_ne _) (hf.ne h₀.prev (G.twin h₀)),
        flipPerm_of_ne (Halfedge.prev_ne _) (Halfedge.prev_ne_next _)
          (hf.ne h₀.prev (G.twin h₀)) (hf.ne h₀.prev (G.twin h₀).next), s1,
        psi_f0 (h₀ := h₀) h₀.prev rfl]
      exact Prod.ext rfl fx.2.2.symm
  by_cases x1 : x.1 = (G.twin h₀).1
  · rcases Halfedge.eq_or_eq_next_or_eq_prev x1 with rfl | rfl | rfl
    · rw [Equiv.swap_apply_right, flipPerm_twin hf, s2, psi_f0 h₀ rfl]
      exact Prod.ext rfl fx.1.symm
    · rw [Equiv.swap_apply_of_ne_of_ne (hf.ne' (G.twin h₀).next h₀) (Halfedge.next_ne _),
        flipPerm_twin_next hf, psi_f1 hf (G.twin h₀).next rfl,
        s_fix _ (Halfedge.ne_next _) (Halfedge.ne_prev _) (hf.ne h₀ (G.twin h₀).next)
          (hf.ne h₀ (G.twin h₀).prev)]
      exact Prod.ext rfl fy.2.1.symm
    · rw [Equiv.swap_apply_of_ne_of_ne (hf.ne' (G.twin h₀).prev h₀) (Halfedge.prev_ne _),
        flipPerm_of_ne (hf.ne' (G.twin h₀).prev h₀) (hf.ne' (G.twin h₀).prev h₀.next)
          (Halfedge.prev_ne _) (Halfedge.prev_ne_next _), s3,
        psi_f1 hf (G.twin h₀).prev rfl]
      exact Prod.ext rfl fy.2.2.symm
  have n0 : x ≠ h₀ := fun e => x0 (by rw [e])
  have n1 : x ≠ G.twin h₀ := fun e => x1 (by rw [e])
  rw [Equiv.swap_apply_of_ne_of_ne n0 n1,
    flipPerm_of_ne n0 (fun e => x0 (by rw [e]; rfl)) n1 (fun e => x1 (by rw [e]; rfl)),
    s_fix x (fun e => x0 (by rw [e]; rfl)) (fun e => x0 (by rw [e]; rfl))
      (fun e => x1 (by rw [e]; rfl)) (fun e => x1 (by rw [e]; rfl)), psi_other x0 x1]

/-- **Flipping the same edge twice gives back the original** up to the explicit isomorphism
that swaps the two faces of the edge and rotates them (`flipFlipShift`). The second flip is at
the new halfedge `next h₀`, which is what the Rust edge id points to after the first flip.

[formula: it-complex-flip] -/
def flipFlipIso : Iso G ((G.flip h₀).flip h₀.next) where
  faces := Equiv.swap h₀.1 (G.twin h₀).1
  shift := G.flipFlipShift h₀
  twin_comm x := by
    change ((G.flip h₀).flip h₀.next).twin (psi G h₀ x) = psi G h₀ (G.twin x)
    have e1 := flipPerm_flipPerm hf (Equiv.swap h₀ (G.twin h₀) x)
    rw [Equiv.swap_apply_self] at e1
    rw [← e1, flip_twin_apply, flip_twin_apply, flipPerm_flipPerm hf, swap_twin_swap]

/-- Hence flipping twice keeps `V`, `E` and `χ` (also immediate from one flip). -/
theorem flip_flip_eulerChar [Fintype F] :
    ((G.flip h₀).flip h₀.next).eulerChar = G.eulerChar :=
  (flipFlipIso hf).eulerChar_eq.symm

end Twice

end Gluing

end DDG.Delaunay
