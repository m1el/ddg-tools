# Integer coordinates in Lean: the face model, and how tier C extends it

Contract: `crates/ddg-integer-coords/SPEC.md` §9. This directory proves tiers A, B and B′; tier C
(the complex model, work package WP9) is meant to be built on the definitions below without
changing them.

## Files

| file | content | ids |
|---|---|---|
| `Formulas.lean` | arithmetic of the oracle's integer formulas (`Ddg.IntegerCoords.*`): signs, Eq. 1–2 on a face given by its arcs (`FaceParams`, `coord_eval`), validity, crossings, Eq. 3 term by term (`flipDoubled_eval`), Eq. 6 | `ic-normal-sign`, `ic-emanating`, `ic-corner-coordinates`, `ic-face-validity`, `ic-reverse-crossing`, `ic-split-edge-normal`, `ic-flip-normal` (identity) |
| `FaceModel.lean` | the face model, its validity, the structure theorem, counts, decoding, the oracle's pairing, Alg. 1 | `ic-face-decoding`, `ic-corner-coordinates`, `ic-emanating`, `ic-cevians-cross` (combinatorial), `ic-trace-step` |
| `FaceSplit.lean` | gaps, separation, the arc families, the oracle's face split (I15) | `ic-face-split-normal` |
| `Flip.lean` | the quadrilateral model `Quad`, Eq. 3 as a count, Eq. 3 applied twice | `ic-flip-normal`, `ic-flip-involution` |
| `VertexRemoval.lean` | Thm. D.1 (angle argument), the failure of the strict form | `ic-vertex-removal` |
| `Chords.lean` | chords of a convex quadrilateral / of a triangle cross iff their ends interleave; cevians | `ic-chord-interleave`, `ic-cevians-cross` |
| `Recover.lean` | App. A (Eq. 15–16), Eq. 5 | `ic-recover-barycentric`, `ic-tangent-length` |
| `Transfer.lean` | Eq. 13–14, the finite part of `P₀ᵀ M_S P₀ = M₀` | `ic-transfer-normal-equations`, `ic-transfer-mass-identity` |
| `SplitEdge.lean` | splitting a crossed edge: `split_side`'s `n_mc` as a separation count | `ic-split-crossed-edge` |
| `SplitAtCrossing.lean` | splitting an edge at its crossing `j`: the arcs crossing the new chord `mc` (`chordCount`, ends interleaving) are `sepCount` at the gap where the arc through `m` does not separate; the Rust `split_side_at_crossing` (`splitSideAtCrossingR`, any side) | `ic-split-at-crossing` |
| `Strip.lean` | `Chord::hit`: the two parameters name one point, uniqueness, both directions | `ic-strip-intersection` |
| `Complex/Basic.lean` | tier C: the complex model `Cx` (gluing, `n` per halfedge, face models = canonical models of the counts), crossings, `rev`, Alg. 1's step `next` = `traceStep` | `ic-trace-step` (on the complex) |
| `Complex/Trace.lean` | walks: end at a corner or return to the start, no crossing twice or in both directions, `2m ≤ Σ n⁺` | `ic-trace-terminates` |
| `Complex/Labels.lean` | corners around a vertex (`rot`), the ends of a corner (`emList`), labels and `RoundaboutsOK` (certify's check), Alg. 3, Alg. 4 | `ic-identify-edge` |
| `Complex/Quad*.lean` | two face models glued on `ij` and re-cut along `kl` (`Glue`): Eq. 3 = arcs crossing `kl` − arcs `k`–`l`, new faces valid, the curves through the quadrilateral unchanged (one file per outer side), the arcs at the corners | `ic-flip-normal` (C part) |
| `Complex/Flip.lean` | the flip of a complex model (`Complex::flip_edge`'s slots) is valid; consistent roundabouts stay consistent with Eq. 4 | `ic-flip-normal`, `ic-roundabout-update` |
| `Complex/Merge.lean` | one side of a straight merge (two faces glued along the spoke) gives a valid face; the Rust `merge_edge_normal` | `ic-merge-edge-normal` |
| `Complex/SplitAtCrossing.lean` | the two faces of a split at a crossing (explicit admissible corner counts), straight with one curve end at `m` (`Merge.Side`, `Merge.Straight`, `Merge.ends = 1`), merging back to the original face; both sides and `Merge.merge_normal` | `ic-split-at-crossing` |
| `Complex/Extract.lean` | every crossing on the curve of exactly one end (walk back), ends are arcs of corners | `ic-extract-intrinsic-edge` (combinatorial) |

Cited (`CITED.md`, "Integer Coordinates"): `ic-normal-isotopy`, `ic-refinement-termination`,
`ic-dijkstra-ball`, `ic-mollify`, `ic-transfer-mass-integral`, `ic-vertex-removal-delta` (open).

## The face model (`FaceModel.lean`)

* A face has corners `0, 1, 2` and sides `0, 1, 2` (`Fin 3`); **side `k` runs from corner `k` to
  corner `k + 1`** (the halfedge `3f + k` of the Rust code) and carries `n k` crossings numbered
  `0 … n k − 1` **from corner `k`** — exactly the paper's combinatorial crossing `(ij→, p)` with `i`
  the tail. Corner `k` lies between side `k + 2` (ending at it) and side `k`; the side opposite
  corner `k` is `k + 1`. `n k` is `n⁺` of the edge (an edge with `n = −1` has no crossings).
* Boundary points `Pt`: `corner k` or `side k p`. `Pt.rank` linearizes the boundary cycle
  `corner 0, side 0, corner 1, side 1, corner 2, side 2` as `0 … N − 1`, `N = n₀ + n₁ + n₂ + 3`.
* A model is `n` and `partner : Fin 3 → ℕ → Pt`, the other end of the arc through each crossing.
  `Valid` = *normal* (side `k` pairs with side `k + 2`, side `k + 1`, or corner `k + 2`), *symmetric*
  (partner of partner), *disjoint* (`Interleave` of the ranks is false for every two arcs; arcs with a
  common corner never interleave). Disjointness is purely combinatorial; `Chords.lean` is the
  geometric justification for straight arcs.
* **Structure theorem** `FaceModel.Valid.eq_canon`: a valid model equals `canon n c` (the oracle's
  `partner`, Alg. 9 / Fig. 11) for admissible corner counts `c` (`Admissible`), and `c` is unique and
  equals Eq. 1 (`Admissible.c_eq`, `Valid.canon_unique`). Everything else about one face is a
  corollary: counts (`Valid.counts`), decoding (`Valid.partner_eq`, `exists_valid_iff`), the oracle
  (`Valid.oracle_partner`), Alg. 1 (`Valid.traceStep_eq`), the arc families (`Valid.arcs`).
* **Regions** (`FaceSplit.lean`) are represented by **gaps** (the boundary piece between ranks `g`
  and `g + 1`); every region of a chord arrangement meets the boundary. `Separates` says an arc
  separates a gap from a boundary point; `sepCount` counts the arcs separating a gap from a corner,
  over the families `cornerArc`/`fanArc`, which enumerate the arcs once each.

## Tier C as built (WP9)

The complex model follows the plan below with one change: the face model of a face is not data but
the canonical model of its counts (`model`, valid by `model_valid`; by `Valid.partner_eq` every
valid model with these counts has the same arcs). Curves are the orbits of `Cx.next`; vertices are
not a type: the labels of an input vertex are functions of its corners that `rot` preserves
(`Labels`, `RoundaboutsOK`). Not modelled: the walk back past inserted vertices (Alg. 3 lines
1–5), the face split as a complex operation, and the flip's effect on the curves at the level of
`Cx.next` (it is stated on the glued quadrilateral, `Glue.Hyp.outer`).

**Build time.** The case analyses are `omega` problems; they are fast only when the facts are linear
and contain no implication: split by the emanating corner of each face first (16 cases), keep
`min`/`max` definitions out of the context unless needed, and try `omega` before `contradiction`.
`Complex/QuadOuter*.lean` are one file per side so that lake builds them in parallel (about a
minute each); the whole IntegerCoords tree builds in about 2.5 minutes from scratch on 64 cores.

## How tier C was planned

1. **Complex model.** A structure with a finite face type `F`, the side gluing
   `twin : F × Fin 3 → F × Fin 3` (an involution without fixed points on interior sides; boundary
   sides unglued), a face model `M f` per face, and the matching conditions
   `(M f).n k = (M f').n k'` for `twin (f, k) = (f', k')`. Crossing `p` of `(f, k)` **is** crossing
   `n − 1 − p` of `(f', k')` (`ic-reverse-crossing`: `Ddg.IntegerCoords.reverseCrossing`). Vertices are
   classes of corners under the gluing; keep the Rust halfedge numbering (`h = 3f + k`) so that the
   oracle's `tracecurve` task and the Rust data structure can be compared literally.
2. **Curves.** A curve step from crossing `(f, k, p)`: pass to the twin `(f', k', n − 1 − p)`, take
   `(M f').partner k' (n − 1 − p)`; a side point continues, a corner ends the curve. This is exactly
   Alg. 1 (`Valid.traceStep_eq` gives the step in the oracle's terms). `ic-trace-terminates`: the
   step map is injective on the finite set of crossings (both `partner` and `twin` are involutions),
   so a trace never repeats a crossing and ends at a corner after at most `Σ n⁺` steps unless it
   closes up; closed curves are excluded by the curves being edges of `T⁰` (they end at vertices), so
   state it for traces starting next to a vertex, or prove "no closed curve" as an invariant.
3. **Flip (`ic-flip-normal` correctness).** Glue the two face models on `ij`; `Quad.ofFaces`
   (`Flip.lean`) already classifies the arcs of the quadrilateral by kind from the two faces' corner
   data (overlaps of the blocks of Fig. 11 along `ij`) and proves the non-crossing constraints
   (`Quad.KL` for both diagonals). Re-cutting along `kl` gives the two new face models whose corner
   data are `Q.flip.faceR`, `Q.flip.faceL` (already proved admissible), and `flip_count` is Eq. 3 as the
   number of arcs crossing `kl`. What remains is to show that the arcs of the glued pair *are* the
   `Quad` kinds (a statement about `partner` maps, by the structure theorem applied to both faces)
   and that the re-cut models glue back to the same curves.
4. **Roundabouts and Alg. 3 (`ic-roundabout-update`, `ic-identify-edge`).** Add, per vertex of `V⁰`,
   the cyclic order of the `T⁰` halfedges and of the `T¹` halfedges around it (a corner cycle as in
   `VertexRemoval.FlatStar`, combinatorial only), with `r` as a map; emanating arcs at a corner are the
   middle block of the opposite side in increasing position (`canon`), which is the ccw order used by
   Eq. 4 and Alg. 3 line 9.
5. **Face split as a complex operation.** `faceSplit?_gapRegion` gives the new counts and
   `splitCount_faceValid` their validity; the remaining step is the gluing of the three new face
   models along the new edges (the arcs of the old face cut by the segments from the region to the
   corners).

## Proof engineering notes

* Rank arithmetic: case on corners with `fin3_cases` (literal `0, 1, 2`) or the `rank_omega`
  macro, rewrite `Fin 3` additions with `fin3_lits`, then `omega`.
* `omega` slows down sharply (or times out in `whnf`) when the context holds many disjunctions or
  hypotheses whose atoms are unfolded definitions; prove small lemmas with minimal contexts
  (`coord_eval`, `FaceParams.rot`), `generalize`/`clear` before `omega`, and use `linarith` for
  linear side conditions in big contexts. `omega` does not split `min` once it has turned it into an
  `if`; rewrite with `min_def` and `split_ifs` first. `omega` does not use a hypothesis `False`;
  `try contradiction` first.
* Build cost: `FaceSplit.lean` takes about 2.5 min (three side lemmas and the validity lemma are
  large case splits), `Flip.lean` about 35 s, the rest seconds.
