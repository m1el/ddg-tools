import DDG.NormalCoords.Complex.Basic
import Mathlib.Tactic.Ring

/-!
# Integer coordinates, tier C: the labels at the input vertices and Alg. 3 (`ic-identify-edge`)

`crates/ddg-integer-coords/SPEC.md` §2.2 (local index, roundabouts, I3), §2.4 item 4 (the
roundabout check of `certify`), §3.1 (Alg. 3, Alg. 4); the Rust `input_directions`,
`identify_at`, `formulas::identify_index`.

**The ends at a corner.** Around a vertex the corners follow each other counter-clockwise by
`rot c = twin (f, k + 2)`: corner `k` of face `f` lies between side `k` (leaving the vertex) and
side `k + 2` (arriving), and the next corner counter-clockwise is the one of the twin of side
`k + 2`, which leaves the vertex (`Complex::outgoing_halfedges`). Inside corner `c = (f, k)`, in
counter-clockwise order, the curves of `T⁰` leave the vertex as follows (`input_directions`): first
side `c` itself if it runs along an input edge (`n = −1`, `nm c = 1`), then the arcs from corner
`k` to the opposite side `k + 1`, in increasing position on that side (`emList c`; in the face
model they are the middle block of side `k + 1`, positions `c_{k+1} … c_{k+1} + e_k − 1`,
`emList_eq`, the Rust `first_emanating`). So corner `c` holds `w c = n⁻ + e` input directions.

**Labels.** The input halfedges at an input vertex `a` are numbered `0 … deg₀(a) − 1`
counter-clockwise (SPEC §2.2, I3), and `r` of a halfedge leaving `a` is the local index of the first
input halfedge at or after it. Since the curves are the input edges drawn on the surface, the
counter-clockwise order of the curve ends at `a` is that of the input halfedges, so the `j`-th
direction of corner `c` has the local index `lab c j = (r c + j) mod deg₀` (no reduction at a
boundary vertex, whose numbering does not wrap: modulus `0`, the oracle's `DEG = 0`), and the
roundabouts of consecutive corners agree with this exactly when the next one is `r + n⁻ + e` (the
check of `certify`, `RoundaboutsOK` below, with the boundary conditions of I3). `RoundaboutsOK` is
the model's form of SPEC §2.4 item 4, first bullet: given it, the labels of all ends around a vertex
are consecutive, so they are the local indices of the input halfedges as soon as one is.

**Results.**

* `Valid.identify_arc` (**`ic-identify-edge`**, Alg. 3 line 9): for the `q`-th arc leaving the
  vertex in corner `c` (first crossing `(f, k + 1, p)`), the oracle's `identifyIndex` of
  `r_ai, n_ai, n_ij, n_ja, p` with `ai = Next(Next(ij))` (the Rust `identify_at`) is its label
  `lab c (n⁻ + q)`; `Valid.identify_along`: an end along side `c` has the label `r c`.
* `Valid.alg4_index` (Alg. 4): the crossing of the `q`-th arc is `p = k + n_ai − r_ai` for its local
  index `k = r_ai + n⁻_ai + q`, as the paper prints it (SPEC §3.1).
* `RoundaboutsOK.lab_next`: the labels continue across corners (`lab` of the first direction
  after corner `c` is `lab c (w c)`), i.e. the labels of the ends are consecutive.

The walk back past inserted vertices on the input edge (Alg. 3 lines 1–5, the Rust loop over
`curve_ends`) is not modelled here; it follows curves through vertices of `V*`, whose two ends are
given by the realization.
-/

namespace DDG.NormalCoords

open DDG.NormalCoords.Formulas (cornerCoord emanating nPlus nMinus identifyIndex)

namespace Cx

variable {F : Type} (C : Cx F)

/-- The next corner counter-clockwise around the vertex of corner `c` (`none` at the boundary). -/
def rot (c : F × Fin 3) : Option (F × Fin 3) := C.twin (c.1, c.2 + 2)

/-- `n⁻` of a side: `1` if it runs along an input edge. -/
def nm (s : F × Fin 3) : ℕ := if C.n s = -1 then 1 else 0

/-- The positions, on the opposite side `k + 1`, of the arcs leaving corner `c = (f, k)`, in
increasing order (counter-clockwise around the corner). -/
def emList (c : F × Fin 3) : List ℕ :=
  (List.range (C.np (c.1, c.2 + 1))).filter
    fun p => decide ((C.M c.1).partner (c.2 + 1) p = .corner c.2)

/-- The number of arcs leaving corner `c`. -/
def emc (c : F × Fin 3) : ℕ := (C.emList c).length

/-- The number of input directions in corner `c`: `n⁻ + e`. -/
def w (c : F × Fin 3) : ℕ := C.nm c + C.emc c

/-- The labels at the input vertices (data of `T⁰` and of the state): which corners are at input
vertices (`V⁰`), `deg₀` and the modulus (`deg₀` at an interior vertex, `0` at a boundary vertex:
`local.modulus`), and the stored roundabout `r` of each halfedge leaving an input vertex. All of
it depends only on the vertex of the corner (`*_rot`). -/
structure Labels (F : Type) where
  input : F × Fin 3 → Prop
  deg : F × Fin 3 → ℕ
  md : F × Fin 3 → ℕ
  r : F × Fin 3 → ℕ

/-- The local index of the `j`-th input direction of corner `c` (`0 ≤ j < w c`). -/
def lab (L : Labels F) (c : F × Fin 3) (j : ℕ) : ℕ := (L.r c + j) % L.md c

/-- **The roundabouts are consistent** (SPEC §2.4 item 4, first bullet; `input_directions`):
labels depend only on the vertex; `r < deg₀`; the next corner's roundabout is `r + n⁻ + e`
(reduced at an interior vertex); at a boundary vertex (modulus `0`, I3) the boundary halfedge
leaving it has `r = 0` and the last corner ends at `deg₀ − 1` (the arriving boundary edge). -/
structure RoundaboutsOK (L : Labels F) : Prop where
  input_rot : ∀ c c', C.rot c = some c' → (L.input c' ↔ L.input c)
  deg_rot : ∀ c c', C.rot c = some c' → L.deg c' = L.deg c
  md_rot : ∀ c c', C.rot c = some c' → L.md c' = L.md c
  md_eq : ∀ c, L.input c → L.md c = 0 ∨ L.md c = L.deg c
  md_bdry : ∀ c, L.input c → C.rot c = none → L.md c = 0
  r_lt : ∀ c, L.input c → L.r c < L.deg c
  r_next : ∀ c c', L.input c → C.rot c = some c' → (L.r c + C.w c) % L.md c = L.r c'
  r_last : ∀ c, L.input c → C.rot c = none → L.r c + C.w c + 1 = L.deg c
  r_first : ∀ c, L.input c → C.twin c = none → L.r c = 0

variable {C}

/-! ## The arcs leaving a corner -/

lemma filter_range_Ico (N a b : ℕ) (hab : a ≤ b) (hb : b ≤ N) :
    (List.range N).filter (fun p => decide (a ≤ p ∧ p < b)) = List.range' a (b - a) := by
  have e : List.range N = List.range' 0 a ++ List.range' a (b - a) ++ List.range' b (N - b) := by
    have h1 := List.range'_append_1 (s := 0) (m := a) (n := b - a)
    have h2 := List.range'_append_1 (s := 0) (m := b) (n := N - b)
    simp only [Nat.zero_add] at h1 h2
    rw [Nat.add_sub_cancel' hab] at h1
    rw [Nat.add_sub_cancel' hb] at h2
    rw [List.range_eq_range', h1, h2]
  rw [e, List.filter_append, List.filter_append, List.filter_eq_nil_iff.2, List.filter_eq_self.2,
    List.filter_eq_nil_iff.2]
  · simp
  · intro x hx; simp only [List.mem_range'_1] at hx; simp only [decide_eq_true_eq]; omega
  · intro x hx; simp only [List.mem_range'_1] at hx; simp only [decide_eq_true_eq]; omega
  · intro x hx; simp only [List.mem_range'_1] at hx; simp only [decide_eq_true_eq]; omega

lemma Valid.admissible (hC : C.Valid) (f : F) : Admissible (C.fn f) (cc (C.fn f)) := by
  apply admissible_cc
  have h := hC.face f
  rw [faceValid_nPlus] at h
  have e : ∀ s, ((C.np s : ℕ) : ℤ) = nPlus (C.n s) := fun s => by simp only [np, nPlus]; omega
  simpa only [fn, e] using h

/-- **The arcs leaving a corner are the middle block of the opposite side**: positions
`c_{k+1}, …, c_{k+1} + e_k − 1` (Fig. 11; the Rust `first_emanating`), with `c` = Eq. 1. -/
lemma Valid.emList_eq (hC : C.Valid) (c : F × Fin 3) :
    C.emList c = List.range' (cc (C.fn c.1) (c.2 + 1))
      (NormalCoords.em (C.fn c.1) (cc (C.fn c.1)) c.2) := by
  have ha := hC.admissible c.1
  have h1 := ha.le (c.2 + 1)
  rw [fin3_add_one_add_one] at h1
  unfold emList M model
  simp only
  rw [List.filter_congr (q := fun p => decide (cc (C.fn c.1) (c.2 + 1) ≤ p ∧
      p < C.fn c.1 (c.2 + 1) - cc (C.fn c.1) (c.2 + 2)))]
  · rw [filter_range_Ico _ _ _ (by omega) (by simp only [fn]; omega)]
    congr 1
    unfold NormalCoords.em; omega
  · intro p hp
    simp only [List.mem_range] at hp
    rw [decide_eq_decide, FaceModel.canon_corner]
    omega

lemma Valid.emc_eq (hC : C.Valid) (c : F × Fin 3) :
    C.emc c = NormalCoords.em (C.fn c.1) (cc (C.fn c.1)) c.2 := by
  rw [emc, hC.emList_eq, List.length_range']

lemma Valid.emList_get (hC : C.Valid) (c : F × Fin 3) {q : ℕ} (hq : q < C.emc c) :
    (C.emList c)[q]'hq = cc (C.fn c.1) (c.2 + 1) + q := by
  simp only [hC.emList_eq c, List.getElem_range'_1]

lemma emanating_nPlus (a b c : ℤ) : emanating a b c = emanating (nPlus a) (nPlus b) (nPlus c) := by
  have h : ∀ x : ℤ, nPlus (nPlus x) = nPlus x := fun x => by simp only [nPlus]; omega
  simp only [emanating, h]

lemma cornerCoord_nPlus (a b c : ℤ) :
    cornerCoord a b c = cornerCoord (nPlus a) (nPlus b) (nPlus c) := by
  have h : ∀ x : ℤ, nPlus (nPlus x) = nPlus x := fun x => by simp only [nPlus]; omega
  simp only [cornerCoord, DDG.NormalCoords.Formulas.cornerNumerator, emanating, h]

lemma np_cast (s : F × Fin 3) : ((C.np s : ℕ) : ℤ) = nPlus (C.n s) := by
  simp only [np, nPlus]; omega

/-- **The number of arcs leaving a corner is Eq. 2** of the face's normal coordinates (what
`input_directions` reads as `coords.e[slot]`). [formula: ic-emanating] -/
theorem Valid.emc_eq_emanating (hC : C.Valid) (c : F × Fin 3) :
    (C.emc c : ℤ) = emanating (C.n (c.1, c.2 + 1)) (C.n c) (C.n (c.1, c.2 + 2)) := by
  rw [hC.emc_eq, ← ((hC.admissible c.1).coords c.2).2, emanating_nPlus (C.n (c.1, c.2 + 1))]
  simp only [fn, np_cast, Prod.mk.eta]

lemma cc_cast (hC : C.Valid) (f : F) (k : Fin 3) :
    (cc (C.fn f) k : ℤ) = cornerCoord (C.n (f, k + 1)) (C.n (f, k)) (C.n (f, k + 2)) := by
  rw [(hC.admissible f).c_eq k, cornerCoord_nPlus (C.n (f, k + 1))]
  simp only [fn, np_cast]

lemma nm_eq (hC : C.Valid) (s : F × Fin 3) : (C.nm s : ℤ) = nMinus (C.n s) := by
  have := hC.n_ge s
  simp only [nm, nMinus]; split_ifs <;> omega

/-- **`ic-identify-edge`: Alg. 3's index is the label of the arc.** The `q`-th arc leaving the
vertex in corner `c = (f, k)` (counter-clockwise; its first crossing is `(ij→, p)` with `ij` side
`k + 1`, `p = (emList c)[q]`) has the label `lab c (n⁻_c + q)`, and that is what Alg. 3 line 9
computes: `identifyIndex r_ai n_ai n_ij n_ja p` with `ai = Next(Next(ij))` = side `k` (the Rust
`identify_at`: `r_ai + n⁻_ai + (p − c^{ja}_i)`, reduced by the modulus). Given `RoundaboutsOK` the
labels are the local indices of the input halfedges (`lab_next`). [formula: ic-identify-edge] -/
theorem Valid.identify_arc (hC : C.Valid) (L : Labels F) (c : F × Fin 3) {q : ℕ}
    (hq : q < C.emc c) :
    identifyIndex (L.r c) (C.n c) (C.n (c.1, c.2 + 1)) (C.n (c.1, c.2 + 2)) ((C.emList c)[q]'hq)
        (L.md c) = (lab L c (C.nm c + q) : ℤ) := by
  rw [hC.emList_get c hq]
  have hcc := cc_cast hC c.1 (c.2 + 1)
  rw [fin3_add_one_add_one, fin3_add_one_add_two] at hcc
  have hnm := nm_eq hC c
  simp only [Prod.mk.eta] at hcc
  unfold identifyIndex lab
  simp only [beq_iff_eq]
  push_cast
  rw [hcc, ← hnm]
  split_ifs with h
  · rw [h, Nat.cast_zero, Int.emod_zero]; ring
  · congr 1; ring

/-- An end along side `c` (`n = −1`) is the first direction of its corner, label `r c` (the Rust
`identify_at`, case `Along`). [formula: ic-identify-edge] -/
theorem identify_along {L : Labels F} (hL : C.RoundaboutsOK L) {c : F × Fin 3} (hc : L.input c) :
    lab L c 0 = L.r c := by
  unfold lab
  rcases hL.md_eq c hc with h | h
  · rw [h, Nat.add_zero, Nat.mod_zero]
  · rw [h, Nat.add_zero, Nat.mod_eq_of_lt (hL.r_lt c hc)]

/-- **Alg. 4's crossing index, as printed** (SPEC §3.1): the `q`-th arc of corner `c` has local
index `k = r_ai + n⁻_ai + q` and first crossing `p = k + n_ai − r_ai`. (It equals
`c^{ja}_i + (k − r_ai − n⁻_ai)` because `c^{ja}_i = n⁺_ai` whenever the corner emanates.) -/
theorem Valid.alg4_index (hC : C.Valid) (L : Labels F) (c : F × Fin 3) {q : ℕ}
    (hq : q < C.emc c) :
    (((C.emList c)[q]'hq : ℕ) : ℤ) = ((L.r c + C.nm c + q : ℕ) : ℤ) + C.n c - L.r c := by
  rw [hC.emList_get c hq]
  have ha := hC.admissible c.1
  have hpos : 0 < NormalCoords.em (C.fn c.1) (cc (C.fn c.1)) c.2 := by
    rw [← hC.emc_eq]; omega
  have h0 := ha.2.1 c.2 hpos
  have h2 := (ha.2.2 c.2 hpos).2
  have hn := ha.n_eq c.2
  have hnm := nm_eq hC c
  have hge := hC.n_ge c
  have hnp := np_cast (C := C) c
  have : cc (C.fn c.1) (c.2 + 1) = C.np c := by
    have e : C.fn c.1 c.2 = C.np c := rfl
    unfold NormalCoords.em at h2
    rw [fin3_add_two_add_one, fin3_add_two_add_two] at h2
    omega
  rw [this]
  simp only [nMinus, nPlus] at hnm hnp
  push_cast
  omega

/-- **The labels are consecutive around a vertex**: the first direction after corner `c` has the
label that a `w c`-th direction of `c` would have. With `lab c (j + 1) = lab c j + 1` (reduced)
inside a corner, the ends around an input vertex get consecutive local indices, which is the model
of the input halfedges' counter-clockwise numbering (SPEC §2.2). [formula: ic-identify-edge] -/
theorem RoundaboutsOK.lab_next {L : Labels F} (hL : C.RoundaboutsOK L) {c c' : F × Fin 3}
    (hc : L.input c) (h : C.rot c = some c') : lab L c' 0 = lab L c (C.w c) := by
  unfold lab
  rw [hL.md_rot c c' h, ← hL.r_next c c' hc h, Nat.add_zero, Nat.mod_mod]

end Cx

end DDG.NormalCoords
