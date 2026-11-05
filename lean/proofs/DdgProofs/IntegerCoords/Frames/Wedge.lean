import DdgProofs.IntegerCoords.Complex.Flip

/-!
# Input-aligned frames: the wedge the roundabout names, and flips (`ic-input-aligned-frame`, part 4 of 4)

`crates/ddg-integer-coords/SPEC.md` §4.12 (I26), §2.2 (local indices, I3), §2.4 item 4; the
Rust `wedge_at` (`crates/ddg-integer-coords/src/geometry.rs`), which `input_strip_of` uses to start
the strip of an edge with `n ≥ 0` at an input vertex, and `frame_directions`
(`crates/ddg-integer-coords/src/frames.rs`). On WP9's complex model (`Cx`, `Labels`,
`RoundaboutsOK`, `Complex/Labels.lean`).

**The wedge.** A halfedge `h` with `n ≥ 0` leaves its input vertex strictly inside an input
corner: no input direction runs along it (`n⁻ = 0`), so the first input direction at or after it
counter-clockwise is an arc leaving the corner, of label `r` (the roundabout; `identify_along`,
`identify_arc`). The input direction just before it is the last direction of the nearest corner
clockwise that holds one (`c₀`, followed by corners without input directions up to the corner of
`h`). `wedge_label`: given the roundabout check (`RoundaboutsOK`), its label is `wedgeIndex r`:
`r − 1` at a boundary vertex (and then `r ≥ 1`), `(r + deg₀ − 1) mod deg₀` inside, and
`r < deg₀`: exactly the local index `wedge_at` computes, and its error branches do not fire. So
`h` leaves into the input corner between the input halfedges of local indices `wedgeIndex r` and
`r`, whose first side is the input halfedge `wedgeIndex r`: the `x` of `Frames/Input.lean`. At a
boundary vertex such a `c₀` always exists: the corner of the boundary halfedge leaving the vertex
runs along the input boundary (`n = −1`) and holds a direction (`w_pos_of_boundary`).

**Flips** (`flip_frame_data`): a flip of another edge keeps, at the moved or kept halfedge
`σ x` and at its twin, the normal coordinate, the gluing, the roundabout and the labels of the tail
(`input`, `deg₀`, the modulus): every integer `frame_directions` reads for an edge that is not the
flipped one. (It also reads the input edge of each crossing, the edge's length and `q⁰` of its
ends, which a flip of another edge does not change either: the lengths and places are not touched,
and the curves through the flipped quadrilateral are the same, `Glue.Hyp.outer`, stated there on
the glued pair, not as `Cx.next` of the flipped complex, which WP9 does not model.)
-/

namespace DDG.IntegerCoords

namespace Cx

variable {F : Type} {C : Cx F}

/-- The local index of the input halfedge before a halfedge of roundabout `r`, as `wedge_at`
computes it: `r − 1` at a boundary vertex (modulus `0`), `(r + deg₀ − 1) mod deg₀` inside. -/
def wedgeIndex (L : Labels F) (c : F × Fin 3) : ℕ :=
  if L.md c = 0 then L.r c - 1 else (L.r c + L.deg c - 1) % L.deg c

/-- Along a chain of corners `cs 0 → cs 1 → … → cs m` counter-clockwise around an input vertex
whose inner corners hold no input direction, the labels stay those of the vertex and the
roundabout is that of the first corner advanced by its directions. -/
theorem RoundaboutsOK.chain {L : Labels F} (hL : C.RoundaboutsOK L) (cs : ℕ → F × Fin 3)
    (m : ℕ) (hrot : ∀ i < m, C.rot (cs i) = some (cs (i + 1)))
    (hw : ∀ i, 0 < i → i < m → C.w (cs i) = 0) (h0 : L.input (cs 0)) :
    ∀ i, 0 < i → i ≤ m → L.input (cs i) ∧ L.md (cs i) = L.md (cs 0) ∧
      L.deg (cs i) = L.deg (cs 0) ∧ L.r (cs i) = (L.r (cs 0) + C.w (cs 0)) % L.md (cs 0) := by
  intro i hi him
  induction i with
  | zero => omega
  | succ k ih =>
    have hr := hrot k (by omega)
    rcases Nat.eq_zero_or_pos k with rfl | hk
    · exact ⟨(hL.input_rot _ _ hr).2 h0, hL.md_rot _ _ hr, hL.deg_rot _ _ hr,
        (hL.r_next _ _ h0 hr).symm⟩
    · obtain ⟨i1, i2, i3, i4⟩ := ih hk (by omega)
      refine ⟨(hL.input_rot _ _ hr).2 i1, by rw [hL.md_rot _ _ hr, i2],
        by rw [hL.deg_rot _ _ hr, i3], ?_⟩
      rw [← hL.r_next _ _ i1 hr, hw k hk (by omega), Nat.add_zero, i4, i2, Nat.mod_mod]

lemma mod_pred (N D : ℕ) (hN : 1 ≤ N) (hD : 1 ≤ D) : (N % D + D - 1) % D = (N - 1) % D := by
  have e := Nat.mod_add_div N D
  rw [← Nat.add_mul_mod_self_left (N % D + D - 1) D (N / D)]
  have : N % D + D - 1 + D * (N / D) = N - 1 + D := by
    generalize N % D = R at e ⊢
    generalize D * (N / D) = P at e ⊢
    omega
  rw [this, Nat.add_mod_right]

/-- **The wedge the roundabout names is the corner `h` leaves into** (`wedge_at`). With
`c = cs m` the corner of `h` (`m ≥ 1`), `c₀ = cs 0` the nearest corner clockwise that holds an
input direction, and no direction in the corners between: the last direction of `c₀`, which is
the input direction just before `h`, has the label `wedgeIndex L c`; at a boundary vertex
`r_c ≥ 1`; and `r_c < deg₀`. [formula: ic-input-aligned-frame] -/
theorem wedge_label {L : Labels F} (hL : C.RoundaboutsOK L) (cs : ℕ → F × Fin 3) {m : ℕ}
    (hm : 0 < m) (hrot : ∀ i < m, C.rot (cs i) = some (cs (i + 1)))
    (hw : ∀ i, 0 < i → i < m → C.w (cs i) = 0) (h0 : L.input (cs 0)) (hw0 : 0 < C.w (cs 0)) :
    lab L (cs 0) (C.w (cs 0) - 1) = wedgeIndex L (cs m) ∧ (L.md (cs m) = 0 → 1 ≤ L.r (cs m)) ∧
      L.r (cs m) < L.deg (cs m) := by
  obtain ⟨i1, i2, i3, i4⟩ := hL.chain cs m hrot hw h0 m hm le_rfl
  have hlt := hL.r_lt _ i1
  refine ⟨?_, ?_, hlt⟩
  · unfold lab wedgeIndex
    rw [i2, i3, i4]
    rcases hL.md_eq _ h0 with hz | hd
    · simp only [hz, ↓reduceIte, Nat.mod_zero]; omega
    · have hD : 1 ≤ L.deg (cs 0) := by have := hL.r_lt _ h0; omega
      have hne : L.md (cs 0) ≠ 0 := by omega
      simp only [hne, ↓reduceIte]
      rw [hd, mod_pred _ _ (by omega) hD]
      congr 1; omega
  · intro hz
    rw [i4, ← i2, hz, Nat.mod_zero]; omega

/-- A halfedge with `n ≥ 0` holds no direction along it (`n⁻ = 0`): the first input direction of
its corner, label `r` (`identify_along`), is an arc leaving the corner (`identify_arc` with
`q = 0`), strictly after `h`; and it is not a boundary halfedge. -/
theorem nm_eq_zero_of_nonneg (hC : C.Valid) {c : F × Fin 3} (hn : 0 ≤ C.n c) :
    C.nm c = 0 ∧ C.twin c ≠ none := by
  refine ⟨by simp only [nm]; split_ifs <;> omega, fun h => ?_⟩
  have := hC.n_bdry c h; omega

/-- At a boundary vertex the corner of the boundary halfedge that leaves it holds an input
direction (the boundary runs along an input edge, `n = −1`): the chain of `wedge_label` starts
there at the latest. -/
theorem w_pos_of_boundary (hC : C.Valid) {c : F × Fin 3} (h : C.twin c = none) : 0 < C.w c := by
  have := hC.n_bdry c h
  simp only [w, nm, this, ↓reduceIte]; omega

namespace FlipAt

variable [DecidableEq F] (D : C.FlipAt)

lemma data_σ (L : Labels F) {x : F × Fin 3} (h1 : x ≠ (D.f, D.s)) (h2 : x ≠ (D.g, D.t)) :
    (C.flip D).n (D.σ x) = C.n x ∧ (C.flip D).twin (D.σ x) = (C.twin x).map D.σ ∧
      ((D.labels L).input (D.σ x) ↔ L.input x) ∧ (D.labels L).deg (D.σ x) = L.deg x ∧
      (D.labels L).md (D.σ x) = L.md x ∧ (D.labels L).r (D.σ x) = L.r x := by
  have a := D.σ_ne_f1 x
  have b := D.σ_ne_g1 x
  refine ⟨?_, ?_, D.labels_σ L h1 h2⟩
  · rw [D.n_of a b, D.σi_σ h1 h2]
  · rw [D.twin_of a b, D.σi_σ h1 h2]

/-- **A flip changes the frame data of no other edge.** For a halfedge `x` not on the flipped
edge, its place `σ x` after the flip (`Complex::flip_edge` moves `j→k`, `i→l` with their data)
has the same normal coordinate, its twin is the place of the old twin, and both keep their
roundabouts and the labels of their tails. These are the integers `frame_directions` reads for
the edge of `x` (`along_direction`, `wedge_at` at both ends); so its directions are unchanged, as
the Rust test checks bitwise (SPEC §4.12, *Invariance*). [formula: ic-input-aligned-frame] -/
theorem flip_frame_data (hC : C.Valid) (L : Labels F) {x : F × Fin 3} (h1 : x ≠ (D.f, D.s))
    (h2 : x ≠ (D.g, D.t)) :
    ((C.flip D).n (D.σ x) = C.n x ∧ (C.flip D).twin (D.σ x) = (C.twin x).map D.σ ∧
      ((D.labels L).input (D.σ x) ↔ L.input x) ∧ (D.labels L).deg (D.σ x) = L.deg x ∧
      (D.labels L).md (D.σ x) = L.md x ∧ (D.labels L).r (D.σ x) = L.r x) ∧
    ∀ z, C.twin x = some z →
      ((D.labels L).input (D.σ z) ↔ L.input z) ∧ (D.labels L).deg (D.σ z) = L.deg z ∧
      (D.labels L).md (D.σ z) = L.md z ∧ (D.labels L).r (D.σ z) = L.r z := by
  refine ⟨D.data_σ L h1 h2, fun z hz => ?_⟩
  obtain ⟨z1, z2⟩ := D.twin_ne_old hC hz h1 h2
  exact D.labels_σ L z1 z2

end FlipAt

end Cx

end DDG.IntegerCoords
