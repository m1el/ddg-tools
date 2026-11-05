import DDG.NormalCoords.Formulas

/-!
# Integer coordinates: the edge flip (Eq. 3) applied twice (tier A)

`crates/ddg-integer-coords/SPEC.md` §4.2, §9.2; the paper Sec. 3.2, Eq. 3.

Edge `ij` with the faces `ijk` and `jil`; the quadrilateral `i, l, j, k` (counter-clockwise). Its arcs
(the pieces of the curves of `T⁰` inside it) are of sixteen kinds (`Quad`): around each corner
(`qi … ql`), between opposite sides (`P1`: `il ↔ jk`, `P2`: `lj ↔ ki`), from a corner to one of the two
sides away from it (`x_i^{lj}`, …), and along a diagonal (`Dij`, an input edge on `ij`, i.e. `n_ij = −1`;
`Dkl`, an arc from `k` to `l`). The two faces `ijk`, `jil` see linear combinations of these
(`Quad.faceR`, `Quad.faceL`); so do the two faces `kil`, `ljk` after the flip, which are the faces of the
relabelled quadrilateral `Quad.flip` (`(i, j, k, l) ↦ (l, k, i, j)`).

* `flip_count` — on a quadrilateral whose arcs do not cross (`Quad.KL`, the pairs of kinds that would
  cross the diagonal `kl` or each other), Eq. 3 of its faces counts the arcs that cross `kl`, minus
  `Dkl` (an arc from `k` to `l` makes the new edge an input edge: `n_kl = −1`).
* `Quad.ofFaces` — two valid faces on `ij` determine the quadrilateral: crossing `p` of `ij` (from `i`)
  is classified by the two faces (Eq. 1–2 blocks of Fig. 11 on either side), and the overlaps of the
  blocks are the counts. The result satisfies `KL` for both diagonals and gives the two old faces back.
* **`ic-flip-involution`** (`flipNormal_flipNormal`) — Eq. 3 applied to the new edge returns `n_ij`,
  for any `n_ij ≥ −1` and any two faces on `ij` that are valid (Eq. 1 numerators even). The other
  edges of the quadrilateral may carry `n = −1` too (Eq. 3 only reads their positive parts).

Notes for the SPEC (reported, not changed here): the paper prints the first term of Eq. 3 as
`c^{jk}_l`, which names no corner (`l` and `jk` share no face); the formula is right with `c^{ji}_l`,
corner `l` of `jil` (as the reference implements, map §2.2), and that is what the oracle computes.
Eq. 3 can return values below `−1` on face pairs that are combinatorially valid but have two or more
arcs from `k` to `l` (impossible for straight arcs in a convex quadrilateral); the involution holds
for them as well.
-/

namespace DDG

namespace NormalCoords

open DDG.NormalCoords.Formulas

/-- The arcs of a quadrilateral `i, l, j, k` by kind. -/
structure Quad where
  qi : ℤ
  qj : ℤ
  qk : ℤ
  ql : ℤ
  /-- between `il` and `jk` -/
  P1 : ℤ
  /-- between `lj` and `ki` -/
  P2 : ℤ
  xilj : ℤ
  xijk : ℤ
  xjki : ℤ
  xjil : ℤ
  xkil : ℤ
  xklj : ℤ
  xljk : ℤ
  xlki : ℤ
  /-- the input edge along `ij` (`n_ij = −1`) -/
  Dij : ℤ
  /-- arcs from `k` to `l` -/
  Dkl : ℤ

namespace Quad

/-- The quadrilateral seen from the other diagonal: `(i, j, k, l) ↦ (l, k, i, j)`. -/
def flip (Q : Quad) : Quad :=
  ⟨Q.ql, Q.qk, Q.qi, Q.qj, Q.P2, Q.P1, Q.xljk, Q.xlki, Q.xkil, Q.xklj, Q.xilj, Q.xijk, Q.xjki,
    Q.xjil, Q.Dkl, Q.Dij⟩

/-- No two arcs cross, as far as the kinds that meet the diagonal `kl` are concerned. -/
def KL (Q : Quad) : Prop :=
  0 ≤ Q.qi ∧ 0 ≤ Q.qj ∧ 0 ≤ Q.qk ∧ 0 ≤ Q.ql ∧ 0 ≤ Q.P1 ∧ 0 ≤ Q.P2 ∧ 0 ≤ Q.xilj ∧ 0 ≤ Q.xijk ∧
    0 ≤ Q.xjki ∧ 0 ≤ Q.xjil ∧ 0 ≤ Q.xkil ∧ 0 ≤ Q.xklj ∧ 0 ≤ Q.xljk ∧ 0 ≤ Q.xlki ∧ 0 ≤ Q.Dij ∧
    0 ≤ Q.Dkl ∧
    (Q.P1 = 0 ∨ Q.P2 = 0) ∧ (Q.P1 = 0 ∨ Q.xklj = 0) ∧ (Q.P1 = 0 ∨ Q.xlki = 0) ∧
    (Q.P2 = 0 ∨ Q.xkil = 0) ∧ (Q.P2 = 0 ∨ Q.xljk = 0) ∧ (Q.xklj = 0 ∨ Q.xljk = 0) ∧
    (Q.xkil = 0 ∨ Q.xlki = 0) ∧
    (Q.Dkl = 0 ∨ (Q.qk = 0 ∧ Q.ql = 0 ∧ Q.P1 = 0 ∧ Q.P2 = 0 ∧ Q.xilj = 0 ∧ Q.xijk = 0 ∧
      Q.xjki = 0 ∧ Q.xjil = 0 ∧ Q.Dij = 0))

/-- The number of arcs crossing the diagonal `kl`, minus `Dkl`. -/
def countKL (Q : Quad) : ℤ :=
  Q.qk + Q.ql + Q.P1 + Q.P2 + Q.xilj + Q.xijk + Q.xjil + Q.xjki + Q.Dij - Q.Dkl

/-- Corner data of the face `ijk`: `(c_i, c_j, c_k, e_i, e_j, e_k)`. -/
def faceR (Q : Quad) : ℤ × ℤ × ℤ × ℤ × ℤ × ℤ :=
  (Q.qi + Q.xlki + Q.P2, Q.qj + Q.xljk + Q.P1, Q.qk, Q.xijk, Q.xjki, Q.xkil + Q.xklj + Q.Dkl)

/-- Corner data of the face `jil`: `(c_j, c_i, c_l, e_j, e_i, e_l)`. -/
def faceL (Q : Quad) : ℤ × ℤ × ℤ × ℤ × ℤ × ℤ :=
  (Q.qj + Q.xklj + Q.P2, Q.qi + Q.xkil + Q.P1, Q.ql, Q.xjil, Q.xilj, Q.xljk + Q.xlki + Q.Dkl)

end Quad

/-- **Eq. 3 counts the arcs crossing the new diagonal.** With the corner data of the two faces of a
quadrilateral whose arcs do not cross, the doubled value of Eq. 3 (as `flipDoubled_eval` writes it)
halved, plus `Dij`, is the number of arcs crossing `kl` minus `Dkl`. [formula: ic-flip-normal] -/
theorem flip_count (Q : Quad) (h : Q.KL) :
    (2 * Q.ql + 2 * Q.qk + iabs ((Q.qj + Q.xklj + Q.P2) - (Q.qj + Q.xljk + Q.P1)) +
        iabs ((Q.qi + Q.xkil + Q.P1) - (Q.qi + Q.xlki + Q.P2)) - (Q.xljk + Q.xlki + Q.Dkl) -
        (Q.xkil + Q.xklj + Q.Dkl) + 2 * (Q.xilj + Q.xijk + Q.xjil + Q.xjki)) / 2 + Q.Dij =
      Q.countKL := by
  obtain ⟨qi, qj, qk, ql, P1, P2, xilj, xijk, xjki, xjil, xkil, xklj, xljk, xlki, Dij, Dkl⟩ := Q
  simp only [Quad.KL, Quad.countKL] at h ⊢
  simp only [iabs]
  split_ifs <;> omega

/-- Length of the overlap of `[a₁, b₁)` and `[a₂, b₂)`. -/
def ov (a1 b1 a2 b2 : ℤ) : ℤ := max 0 (min b1 b2 - max a1 a2)

/-- The quadrilateral of two faces on `ij`: crossing `p` of `ij` (from `i`) is, in `ijk`, in the
block `[0, c_i)` (around `i`), `[c_i, c_i + e_k)` (from `k`) or `[c_i + e_k, n)` (around `j`), and in
`jil` in `[0, c'_i)`, `[c'_i, c'_i + e_l)` or `[c'_i + e_l, n)`; the overlaps of the blocks count the
arcs through `ij` by kind. The other kinds are corner data of one face. -/
def Quad.ofFaces (ci cj ck ei ej ek ci' cl ej' ei' el m : ℤ) : Quad where
  qi := ov 0 ci 0 ci'
  qj := ov (ci + ek) (ci + cj + ek) (ci' + el) (ci + cj + ek)
  qk := ck
  ql := cl
  P1 := ov (ci + ek) (ci + cj + ek) 0 ci'
  P2 := ov 0 ci (ci' + el) (ci + cj + ek)
  xilj := ei'
  xijk := ei
  xjki := ej
  xjil := ej'
  xkil := ov ci (ci + ek) 0 ci'
  xklj := ov ci (ci + ek) (ci' + el) (ci + cj + ek)
  xljk := ov (ci + ek) (ci + cj + ek) ci' (ci' + el)
  xlki := ov 0 ci ci' (ci' + el)
  Dij := m
  Dkl := ov ci (ci + ek) ci' (ci' + el)

set_option maxHeartbeats 4000000 in
/-- **Two valid faces on `ij` make a quadrilateral without crossings**, whose faces are the given
ones and whose faces after the flip are valid. [formula: ic-flip-involution] -/
theorem Quad.ofFaces_spec {ci cj ck ei ej ek cj' ci' cl ej' ei' el m : ℤ}
    (p1 : FaceParams ci cj ck ei ej ek) (p2 : FaceParams cj' ci' cl ej' ei' el)
    (hA : ci + cj + ek = cj' + ci' + el) (hm0 : 0 ≤ m) (hm1 : m ≤ 1)
    (hmA : m = 1 → ci + cj + ek = 0) :
    let Q := Quad.ofFaces ci cj ck ei ej ek ci' cl ej' ei' el m
    Q.faceR = (ci, cj, ck, ei, ej, ek) ∧ Q.faceL = (cj', ci', cl, ej', ei', el) ∧ Q.Dij = m ∧
      Q.KL ∧ Q.flip.KL ∧
      FaceParams Q.flip.faceR.1 Q.flip.faceR.2.1 Q.flip.faceR.2.2.1 Q.flip.faceR.2.2.2.1
        Q.flip.faceR.2.2.2.2.1 Q.flip.faceR.2.2.2.2.2 ∧
      FaceParams Q.flip.faceL.1 Q.flip.faceL.2.1 Q.flip.faceL.2.2.1 Q.flip.faceL.2.2.2.1
        Q.flip.faceL.2.2.2.2.1 Q.flip.faceL.2.2.2.2.2 := by
  intro Q
  unfold FaceParams at p1 p2
  simp only [Q, Quad.ofFaces, Quad.faceR, Quad.faceL, Quad.flip, Quad.KL, FaceParams, Prod.mk.injEq]
  refine ⟨⟨?_, ?_, ?_, ?_, ?_, ?_⟩, ⟨?_, ?_, ?_, ?_, ?_, ?_⟩, trivial, ?_, ?_, ?_, ?_⟩
  all_goals first | rfl | (simp only [ov] <;> omega)

/-- The same, as an existence statement about an opaque quadrilateral. -/
theorem Quad.exists_of_faces {ci cj ck ei ej ek cj' ci' cl ej' ei' el m : ℤ}
    (p1 : FaceParams ci cj ck ei ej ek) (p2 : FaceParams cj' ci' cl ej' ei' el)
    (hA : ci + cj + ek = cj' + ci' + el) (hm0 : 0 ≤ m) (hm1 : m ≤ 1)
    (hmA : m = 1 → ci + cj + ek = 0) :
    ∃ Q : Quad, Q.faceR = (ci, cj, ck, ei, ej, ek) ∧ Q.faceL = (cj', ci', cl, ej', ei', el) ∧
      Q.Dij = m ∧ Q.KL ∧ Q.flip.KL ∧
      FaceParams Q.flip.faceR.1 Q.flip.faceR.2.1 Q.flip.faceR.2.2.1 Q.flip.faceR.2.2.2.1
        Q.flip.faceR.2.2.2.2.1 Q.flip.faceR.2.2.2.2.2 ∧
      FaceParams Q.flip.faceL.1 Q.flip.faceL.2.1 Q.flip.faceL.2.2.1 Q.flip.faceL.2.2.2.1
        Q.flip.faceL.2.2.2.2.1 Q.flip.faceL.2.2.2.2.2 :=
  ⟨_, Quad.ofFaces_spec p1 p2 hA hm0 hm1 hmA⟩

/-- **`ic-flip-involution`: Eq. 3 applied twice returns `n_ij`.** For an edge `ij` with
`n_ij ≥ −1` whose two faces `ijk` (`n_ij, n_jk, n_ki`) and `jil` (`n_ij, n_il, n_lj`) are valid
(`faceValid`, i.e. the numerators of Eq. 1 are even: the counts come from a pair of face models),
flipping gives `n_kl = flipNormal n_ij n_jk n_ki n_il n_lj`, and flipping `kl` back, with the faces
`lki` and `klj` (the relabelling `(i, j, k, l) ↦ (l, k, i, j)`), gives `n_ij` again. The other four
normal coordinates may be `−1` (shared edges) or anything else; only their positive parts enter.
[formula: ic-flip-involution] -/
theorem flipNormal_flipNormal {a b c d e : ℤ} (ha : -1 ≤ a) (v1 : faceValid a b c = true)
    (v2 : faceValid a d e = true) : flipNormal (flipNormal a b c d e) c d e b = a := by
  obtain ⟨ci, cj, ck, ei, ej, ek, p1, h_ij, h_jk, h_ki⟩ := params_of_faceValid v1
  obtain ⟨cj', ci', cl, ej', ei', el, p2, h_ij', h_il, h_lj⟩ := params_of_faceValid v2
  have hpa : nPlus a = if a = -1 then 0 else a := by simp only [nPlus]; split_ifs <;> omega
  have hma : nMinus a = if a = -1 then 1 else 0 := by simp only [nMinus]; split_ifs <;> omega
  obtain ⟨Q, hR, hL, hD, hKL, hKL', pR', pL'⟩ := Quad.exists_of_faces (m := nMinus a) p1 p2
    (by rw [← h_ij, ← h_ij']) (by rw [hma]; split_ifs <;> omega) (by rw [hma]; split_ifs <;> omega)
    (by intro h; rw [← h_ij, hpa]; rw [hma] at h; split_ifs at h ⊢ <;> omega)
  simp only [Quad.faceR, Quad.faceL, Prod.mk.injEq] at hR hL
  obtain ⟨rfl, rfl, rfl, rfl, rfl, rfl⟩ := hR
  obtain ⟨rfl, rfl, rfl, rfl, rfl, rfl⟩ := hL
  -- the first flip counts the arcs crossing `kl`
  have hX : flipNormal a b c d e = Q.countKL := by
    rw [flipNormal, flipDoubled_eval p1 p2 h_ij h_jk h_ki h_ij' h_il h_lj, ← hD]
    exact flip_count Q hKL
  rw [hX]
  -- the second flip, on the relabelled quadrilateral
  have hK := hKL
  simp only [Quad.KL] at hK
  simp only [Quad.flip, Quad.faceR, Quad.faceL] at pR' pL'
  have np : nPlus Q.countKL = (Q.ql + Q.xjil + Q.P1) + (Q.qk + Q.xjki + Q.P2) +
      (Q.xilj + Q.xijk + Q.Dij) := by
    simp only [nPlus, Quad.countKL]; omega
  have nm : nMinus Q.countKL = Q.Dkl := by
    simp only [nMinus, Quad.countKL]; omega
  have e2 := flipDoubled_eval pR' pL' np (by rw [h_ki]; linarith) (by rw [h_il]; linarith)
    (by rw [np]; linarith) (by rw [h_lj]; linarith) (by rw [h_jk]; linarith)
  have f2 := flip_count Q.flip hKL'
  simp only [Quad.flip, Quad.countKL] at f2
  rw [flipNormal, e2, nm]
  rw [f2]
  have : a = nPlus a - nMinus a := (nPlus_sub_nMinus a).symm
  rw [this, h_ij, ← hD]
  linarith

end NormalCoords

end DDG
