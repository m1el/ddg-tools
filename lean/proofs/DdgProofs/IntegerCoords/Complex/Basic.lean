import DdgProofs.IntegerCoords.FaceModel
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic.Abel

/-!
# Integer coordinates, tier C: the complex model

`crates/ddg-integer-coords/SPEC.md` §9.1 ("a complex model"), §2.1–2.4; `DESIGN.md` in the parent
directory, "How tier C should extend it". Work package WP9.

A **complex model** is a finite Δ-complex given by gluing data and one normal coordinate per side:

* faces `F` (a finite type), each with corners `0, 1, 2` and sides `0, 1, 2`; side `k` of face `f`
  runs from corner `k` to corner `k + 1` (the halfedge `h = 3f + k` of the Rust code, which is the
  only numbering used below: a side **is** a halfedge, `(f, k) ↔ 3f + k`);
* `twin : F × Fin 3 → Option (F × Fin 3)`, the side pairing (`Complex::twin`; `none` on the
  boundary), an involution without fixed points; gluing reverses the orientation, which is built
  into the indexing (side `(f, k)` from corner `k` to `k + 1` is glued to `(f', k')` from `k'` to
  `k' + 1`, so corner `k` of `f` is corner `k' + 1` of `f'`). Nothing else is required: a face may
  be glued to itself and use a vertex several times, as the Rust `Complex` allows (SPEC §2.1);
* `n : F × Fin 3 → ℤ`, the normal coordinate of the edge of each halfedge (`halfedge_normal`),
  equal on twins, `≥ −1`, and `−1` on boundary edges (SPEC §2.4, "ranges").

`Cx.Valid` also asks that every face be valid in the sense of Eq. 1 (the oracle's `faceValid`;
SPEC §2.4 item 1). The **face model** of a face is then not extra data: by the structure theorem of
tier B (`FaceModel.Valid.eq_canon`, `FaceModel.Valid.partner_eq`) a valid face model is determined
by its three side counts and is the canonical pairing `canon n c` with `c` = Eq. 1. So `Cx.M f` is
that canonical model (`model`), and `model_valid` proves it valid; any other valid model with the
same counts has the same arcs. The matching of crossing counts on glued sides (§9.1) is `n_twin`;
crossing `p` of a side is crossing `n⁺ − 1 − p` of its twin (`ic-reverse-crossing`), which is how
`rev` and `next` below glue the face models.

**Crossings.** A crossing `z = ((f, k), p)` with `p < n⁺` is the Rust `Crossing { halfedge: 3f + k,
index: p }`: the curve leaves face `f` through side `k` at its `p`-th crossing from the tail
(SPEC §3.1, the convention of Supp. Alg. 1). `rev z` is the same point crossed the other way
(`Crossing::reversed`), and `next z` is one step of Alg. 1 (`trace_with`): into the face of the twin,
along the arc of its face model through the entry point, to the next crossing (`inl`) or to a
corner, where the curve ends (`inr`, the Rust `Arrival::corner`). `next_eq_traceStep` shows that
this is exactly what the oracle's (and the Rust) `traceStep` computes from the integers.

The **curves** of §9.1 (maximal chains of arcs glued through sides, ending at corners) are the
orbits of `next`; `DdgProofs/IntegerCoords/Complex/Trace.lean` proves what `TraceFrom` does with
them (`ic-trace-terminates`). Disjointness of arcs is inherited from the face models
(combinatorial, by non-interleaving: no Jordan curve theorem is needed anywhere).
-/

namespace DDG.IntegerCoords

open Ddg.IntegerCoords (cornerCoord cornerNumerator emanating faceValid nPlus nMinus traceStep)

/-! ## The face model of side counts -/

/-- Eq. 1 at every corner, as natural numbers: the corner counts of the canonical model. -/
def cc (n : Fin 3 → ℕ) (k : Fin 3) : ℕ := (cornerCoord (n (k + 1)) (n k) (n (k + 2))).toNat

/-- The face model of side counts `n`: the canonical pairing for the corner counts of Eq. 1
(Fig. 11, Alg. 9; the oracle's `partner`). -/
def model (n : Fin 3 → ℕ) : FaceModel := ⟨n, canon n (cc n)⟩

lemma faceValid_nPlus (a b c : ℤ) : faceValid a b c = faceValid (nPlus a) (nPlus b) (nPlus c) := by
  have h : ∀ x : ℤ, nPlus (nPlus x) = nPlus x := fun x => by simp only [nPlus]; omega
  simp only [faceValid, cornerNumerator, emanating, h]

/-- The corner counts of a valid face are admissible: Eq. 1 is the canonical model's `c`. -/
lemma admissible_cc {n : Fin 3 → ℕ} (h : faceValid (n 0) (n 1) (n 2) = true) :
    Admissible n (cc n) := by
  obtain ⟨M, hMn, hM⟩ := (FaceModel.exists_valid_iff n).2 h
  obtain ⟨c, hc, -⟩ := hM.eq_canon
  rw [hMn] at hc
  have : cc n = c := by
    funext k
    have := hc.c_eq k
    simp only [cc, ← this, Int.toNat_natCast]
  rwa [this]

/-- **The face model of a valid face is valid**; it is the only valid model with these side
counts (`FaceModel.Valid.partner_eq`). [formula: ic-face-decoding] -/
theorem model_valid {n : Fin 3 → ℕ} (h : faceValid (n 0) (n 1) (n 2) = true) : (model n).Valid :=
  canon_valid (admissible_cc h)

/-! ## Complex models -/

/-- A complex model: the side pairing (`none` on the boundary) and the normal coordinate of the
edge of every side (halfedge `3f + k ↔ (f, k)`). -/
structure Cx (F : Type) where
  /-- the twin of a side (`Complex::twin`) -/
  twin : F × Fin 3 → Option (F × Fin 3)
  /-- the normal coordinate of the edge of a side -/
  n : F × Fin 3 → ℤ

/-- A crossing: a side and a position on it, counted from the side's tail. -/
abbrev Xing (F : Type) := (F × Fin 3) × ℕ

namespace Cx

variable {F : Type} (C : Cx F)

/-- `n⁺` of a side, as a natural number. -/
def np (s : F × Fin 3) : ℕ := (C.n s).toNat

/-- The side counts of face `f`. -/
def fn (f : F) : Fin 3 → ℕ := fun k => C.np (f, k)

/-- The face model of face `f`. -/
def M (f : F) : FaceModel := model (C.fn f)

/-- Validity of a complex model: the gluing is an involution without fixed points, the normal
coordinates are `≥ −1`, equal on twins and `−1` on the boundary (SPEC §2.4 "ranges"), and every
face is valid (Eq. 1 numerators even, SPEC §2.4 item 1). -/
structure Valid : Prop where
  twin_symm : ∀ s t, C.twin s = some t → C.twin t = some s
  twin_ne : ∀ s, C.twin s ≠ some s
  n_ge : ∀ s, -1 ≤ C.n s
  n_twin : ∀ s t, C.twin s = some t → C.n t = C.n s
  n_bdry : ∀ s, C.twin s = none → C.n s = -1
  face : ∀ f, faceValid (C.n (f, 0)) (C.n (f, 1)) (C.n (f, 2)) = true

/-- `z` is a crossing: its position is below `n⁺` of its side. -/
def IsX (z : Xing F) : Prop := z.2 < C.np z.1

/-- The twin, or the side itself on the boundary (where it is never used). -/
def twinD (s : F × Fin 3) : F × Fin 3 := (C.twin s).getD s

/-- The same crossing seen from the twin (`Crossing::reversed`, `ic-reverse-crossing`). -/
def rev (z : Xing F) : Xing F := (C.twinD z.1, C.np z.1 - 1 - z.2)

/-- One step of Supp. Alg. 1 (`trace_with`): from crossing `z`, into the face of the twin at the
reversed position, along the arc of its face model: the next crossing (`inl`), or the corner where
the curve ends (`inr`). -/
def next (z : Xing F) : Xing F ⊕ (F × Fin 3) :=
  match (C.M (C.twinD z.1).1).partner (C.twinD z.1).2 (C.np z.1 - 1 - z.2) with
  | .side j r => .inl (((C.twinD z.1).1, j), r)
  | .corner c => .inr ((C.twinD z.1).1, c)

/-- `z` is an **end** of a curve at a vertex: the arc of its own face through it comes from a
corner (the first crossing of a curve leaving a vertex, SPEC §3.1). -/
def IsEnd (z : Xing F) : Prop := ∃ c, (C.M z.1.1).partner z.1.2 z.2 = .corner c

variable {C}

lemma M_n (f : F) (k : Fin 3) : (C.M f).n k = C.np (f, k) := rfl

lemma Valid.M_valid (hC : C.Valid) (f : F) : (C.M f).Valid := by
  apply model_valid
  have h := hC.face f
  rw [faceValid_nPlus] at h
  have e : ∀ s, ((C.np s : ℕ) : ℤ) = nPlus (C.n s) := fun s => by
    simp only [np, nPlus]; omega
  simpa only [fn, e] using h

lemma Valid.twin_of_pos (hC : C.Valid) {s : F × Fin 3} (h : 0 < C.np s) :
    ∃ t, C.twin s = some t := by
  rcases e : C.twin s with _ | t
  · have := hC.n_bdry s e; simp only [np] at h; omega
  · exact ⟨t, rfl⟩

lemma Valid.twinD_eq (hC : C.Valid) {s : F × Fin 3} (h : 0 < C.np s) :
    C.twin s = some (C.twinD s) := by
  obtain ⟨t, ht⟩ := hC.twin_of_pos h
  simp [twinD, ht]

lemma Valid.np_twinD (hC : C.Valid) {s : F × Fin 3} (h : 0 < C.np s) :
    C.np (C.twinD s) = C.np s := by
  simp only [np, hC.n_twin _ _ (hC.twinD_eq h)]

lemma Valid.twinD_twinD (hC : C.Valid) {s : F × Fin 3} (h : 0 < C.np s) :
    C.twinD (C.twinD s) = s := by
  have h2 := hC.twin_symm _ _ (hC.twinD_eq h)
  show (C.twin (C.twinD s)).getD (C.twinD s) = s
  rw [h2]; rfl

lemma Valid.rev_isX (hC : C.Valid) {z : Xing F} (hz : C.IsX z) : C.IsX (C.rev z) := by
  simp only [IsX, rev] at *
  rw [hC.np_twinD (by omega)]; omega

lemma Valid.rev_rev (hC : C.Valid) {z : Xing F} (hz : C.IsX z) : C.rev (C.rev z) = z := by
  simp only [IsX] at hz
  simp only [rev, hC.twinD_twinD (by omega : 0 < C.np z.1), hC.np_twinD (by omega : 0 < C.np z.1)]
  ext <;> simp; omega

lemma Valid.rev_ne (hC : C.Valid) {z : Xing F} (hz : C.IsX z) : C.rev z ≠ z := by
  intro h
  have h1 := hC.twinD_eq (s := z.1) (by simp only [IsX] at hz; omega)
  have : C.twinD z.1 = z.1 := congrArg Prod.fst h
  rw [this] at h1
  exact hC.twin_ne _ h1

/-! ### Steps -/

lemma next_eq_inl {z y : Xing F} (h : C.next z = .inl y) :
    (C.M (C.twinD z.1).1).partner (C.twinD z.1).2 (C.np z.1 - 1 - z.2) = .side y.1.2 y.2 ∧
      y.1.1 = (C.twinD z.1).1 := by
  unfold next at h
  split at h
  · rename_i j r e; cases h; exact ⟨e, rfl⟩
  · cases h

lemma next_eq_inr {z : Xing F} {c : F × Fin 3} (h : C.next z = .inr c) :
    (C.M (C.twinD z.1).1).partner (C.twinD z.1).2 (C.np z.1 - 1 - z.2) = .corner c.2 ∧
      c.1 = (C.twinD z.1).1 := by
  unfold next at h
  split at h
  · cases h
  · rename_i c' e; cases h; exact ⟨e, rfl⟩

/-- The next crossing is a crossing. -/
lemma Valid.next_isX (hC : C.Valid) {z y : Xing F} (hz : C.IsX z) (h : C.next z = .inl y) :
    C.IsX y := by
  obtain ⟨e, hy⟩ := next_eq_inl h
  have hp : C.np z.1 - 1 - z.2 < (C.M (C.twinD z.1).1).n (C.twinD z.1).2 := by
    rw [M_n, Prod.mk.eta, hC.np_twinD (by simp only [IsX] at hz; omega)]
    simp only [IsX] at hz; omega
  have := ((hC.M_valid _).side_lt hp e).1
  rw [M_n, ← hy] at this
  exact this

/-- **The step map commutes with reversal**: if the curve goes from `z` to `y`, then walked
backwards it goes from `rev y` to `rev z`. -/
lemma Valid.next_rev (hC : C.Valid) {z y : Xing F} (hz : C.IsX z) (h : C.next z = .inl y) :
    C.next (C.rev y) = .inl (C.rev z) := by
  have hy := hC.next_isX hz h
  obtain ⟨e, hy1⟩ := next_eq_inl h
  simp only [IsX] at hz hy
  have hpz : 0 < C.np z.1 := by omega
  have hpy : 0 < C.np y.1 := by omega
  set t := C.twinD z.1 with ht
  have hp : C.np z.1 - 1 - z.2 < (C.M t.1).n t.2 := by
    rw [M_n, Prod.mk.eta, ht, hC.np_twinD hpz]; omega
  have hs := (hC.M_valid t.1).symm t.2 _ y.1.2 y.2 hp e
  -- the step from `rev y`
  have hty : C.twinD (C.twinD y.1) = y.1 := hC.twinD_twinD hpy
  have hy' : y.1 = (t.1, y.1.2) := by rw [← hy1]
  unfold next rev
  simp only [hty, hC.np_twinD hpy]
  have hq : C.np y.1 - 1 - (C.np y.1 - 1 - y.2) = y.2 := by omega
  rw [hq]
  conv => lhs; rw [hy']
  simp only
  rw [hs]

/-- A crossing never steps to its own reversal (the partner of a crossing is on another side). -/
lemma Valid.next_ne_rev (hC : C.Valid) {z : Xing F} (hz : C.IsX z) :
    C.next z ≠ .inl (C.rev z) := by
  intro h
  obtain ⟨e, -⟩ := next_eq_inl h
  simp only [IsX] at hz
  have hp : C.np z.1 - 1 - z.2 < (C.M (C.twinD z.1).1).n (C.twinD z.1).2 := by
    rw [M_n, Prod.mk.eta, hC.np_twinD (by omega)]; omega
  have := ((hC.M_valid _).side_lt hp e).2
  simp only [rev] at this
  rcases this with h1 | h1
  · exact fin3_add_two_ne _ h1.symm
  · exact fin3_add_one_ne _ h1.symm

/-- **The step map is injective** on crossings. -/
lemma Valid.next_inj (hC : C.Valid) {z z' y : Xing F} (hz : C.IsX z) (hz' : C.IsX z')
    (h : C.next z = .inl y) (h' : C.next z' = .inl y) : z = z' := by
  have e1 := hC.next_rev hz h
  have e2 := hC.next_rev hz' h'
  rw [e1] at e2
  have e3 : C.rev z = C.rev z' := Sum.inl.inj e2
  rw [← hC.rev_rev hz, e3, hC.rev_rev hz']

/-- **Nothing steps to an end**: the crossing through which an arc leaves a corner is the first
crossing of its curve. -/
lemma Valid.not_next_eq_end (hC : C.Valid) {x z : Xing F} (hx : C.IsX x) (he : C.IsEnd z) : C.next x ≠ .inl z := by
  intro h
  obtain ⟨c, hc⟩ := he
  obtain ⟨e, hz1⟩ := next_eq_inl h
  simp only [IsX] at hx
  have hp : C.np x.1 - 1 - x.2 < (C.M (C.twinD x.1).1).n (C.twinD x.1).2 := by
    rw [M_n, Prod.mk.eta, hC.np_twinD (by omega)]; omega
  have hs := (hC.M_valid _).symm _ _ _ _ hp e
  rw [← hz1, hc] at hs
  cases hs

/-! ### The step computed from the integers (Supp. Alg. 1, the oracle's `traceStep`) -/

lemma traceStep_nPlus (a b c p : ℤ) :
    traceStep a b c p = traceStep (nPlus a) (nPlus b) (nPlus c) p := by
  have h : ∀ x : ℤ, nPlus (nPlus x) = nPlus x := fun x => by simp only [nPlus]; omega
  simp only [traceStep, cornerCoord, cornerNumerator, emanating, h]
  rfl

/-- What the Rust `trace_with` (and the oracle's `traceStep`) does with crossing `z`. -/
def oracleNext (z : Xing F) : Xing F ⊕ (F × Fin 3) :=
  let t := C.twinD z.1
  match traceStep (C.n z.1) (C.n (t.1, t.2 + 1)) (C.n (t.1, t.2 + 2)) z.2 with
  | .toIK q => .inl ((t.1, t.2 + 1), q.toNat)
  | .toKJ q => .inl ((t.1, t.2 + 2), q.toNat)
  | .vertex => .inr (t.1, t.2 + 2)

/-- **The model's step is Alg. 1.** On every crossing of a valid complex model, `next` is the step
that the oracle's `traceStep` computes from the three normal coordinates of the face entered and
the crossing index (the Rust `trace_with` makes the same three cases). [formula: ic-trace-step] -/
theorem Valid.next_eq_oracle (hC : C.Valid) {z : Xing F} (hz : C.IsX z) :
    C.next z = C.oracleNext z := by
  simp only [IsX] at hz
  have hpz : 0 < C.np z.1 := by omega
  set t := C.twinD z.1 with ht
  have hn : ∀ s, ((C.np s : ℕ) : ℤ) = nPlus (C.n s) := fun s => by simp only [np, nPlus]; omega
  have hnt : C.np t = C.np z.1 := hC.np_twinD hpz
  have hst := (hC.M_valid t.1).traceStep_eq t.2 z.2 (by rw [M_n, Prod.mk.eta, hnt]; exact hz)
  simp only [M_n, Prod.mk.eta] at hst
  have hnz : nPlus (C.n z.1) = nPlus (C.n t) := by rw [← hn, ← hn, hnt]
  have key : traceStep (C.n z.1) (C.n (t.1, t.2 + 1)) (C.n (t.1, t.2 + 2)) z.2 =
      traceStep (C.np t) (C.np (t.1, t.2 + 1)) (C.np (t.1, t.2 + 2)) z.2 := by
    rw [traceStep_nPlus, hnz, ← hn, ← hn, ← hn]
  rw [hnt] at hst key
  unfold next oracleNext
  dsimp only
  rw [← ht, key, ← hst]
  rcases traceStep (C.np z.1) (C.np (t.1, t.2 + 1)) (C.np (t.1, t.2 + 2)) z.2 with q | q | - <;>
    simp [stepPt]

end Cx

end DDG.IntegerCoords
