# Delaunay layer 4: discharging the finiteness hypothesis

Goal 1 of layer 4: prove `MetricGluing.BoundedLengthFinite s` for every `s : MetricGluing F`
(`DdgProofs/Delaunay/Termination.lean`), so that intrinsic Delaunay flipping terminates
unconditionally and a Delaunay triangulation with the original vertices exists unconditionally
(Indermitte–Liebling–Troyanov–Clémençon 2001, Thm. 1; Bobenko–Springborn 2007, Prop. 12 /
Thm. 15; Fisher–Springborn–Bobenko–Schröder 2007).

## The route: an explicit overlay of the flipped triangulation on the original one

Three routes were compared.

* (a) *The glued surface as a metric space.* Build `|s|` as a quotient of Euclidean triangles,
  show it is a length space, define geodesics, prove finiteness of bounded geodesics. Nothing of
  this exists in Mathlib (no piecewise flat surfaces, no quotient length metrics suited to this,
  no "locally geodesic" notion without charts); the statement "an edge of a flipped triangulation
  is a geodesic of `|s|`" needs charts anyway. Far too expensive; rejected.
* (b) *Normal coordinates / crossing sequences.* Track every edge of `t` as the sequence of edges
  of `s` it crosses. `DdgProofs/IntegerCoords/` has the combinatorial face models and the
  combinatorial flip rule, but no geometry: that the straight diagonal of the flipped diamond
  crosses the arcs as the rule predicts, and that the length is the developed distance, would
  have to be proved — and that needs to know how the edges of `s` run *inside* each face of `t`,
  i.e. the overlay. So (b) reduces to (c) plus a combinatorial layer that is not needed.
* (c) **Chosen: the overlay as planar cells, with pointwise "germ" consistency.** For each face
  `τ` of `t` keep a finite set of *cells*: faces `σ` of `s` with a rigid motion placing the
  standard layout of `σ` into the plane in which `τ` is laid out. The invariants are local and
  planar: the cells cover `τ`; at every interior point of `τ`, and at every point of an open
  side of `τ` seen together with the neighbouring face placed across that side, the cells
  containing the point form a *consistent germ* (one cell with the point interior, or two cells
  that are the unfolding of each other across an interior edge of `s` with the point on that
  open side; on a boundary edge of `t`, one cell with the point on a boundary side of `s`); the
  corners of `τ` that lie in a cell are corners of that cell. **A flip changes nothing in the
  cells**: the two faces of the diamond are laid out as one convex quadrilateral `D`; the cells of
  the two new faces are the cells of `D` that overlap the new face. Every edge of `t` is then a
  straight segment through a chain of cells of `s`, developed in the plane; its length is a
  function of the combinatorial *word* (first face, first corner, the exit sides crossed, last
  corner); and a segment of length `≤ L` crosses at most `N(s, L)` sides, by a planar argument
  (around one vertex the crossed corner angles add up to less than `π`; switching the pivot
  vertex costs at least `sin α_min · ℓ_min` of length). Finitely many words of bounded length,
  finitely many gluings on the fixed face type: `BoundedLengthFinite` follows (even without the
  "up to isomorphism").

Why (c) is the cheapest honest route: it needs no topology of surfaces, no Jordan curve theorem,
no quotient spaces, no manifold charts, and no new combinatorial theory. Every lemma is about
triangles in `ℝ × ℝ` (`sda`, `dist2`, `rigid`, hull membership by signs — all in
`DdgProofs/Delaunay/PlanarBasic.lean`, `DdgProofs/Intrinsic/Diamond.lean`,
`DdgProofs/Bff/Quality.lean`) plus finite sets. The phantom-placement problem (developing around
a cone point does not return to the same placement) is handled by triggering every consistency
condition at a *point* of the face or of its open sides, never by "the placed triangles overlap".

What is **not** modelled, and why it is not needed: the surface itself, geodesics as curves,
homotopy classes, uniqueness of geodesics. The theorem proved is literally the hypothesis
`BoundedLengthFinite` of layer 3, for every state, closed or with boundary, with self-loops,
multiple edges and self-glued faces allowed.

## Stages (files under `DdgProofs/Delaunay/Surface/`)

| Stage | File | Content | Est. lines | Status |
|---|---|---|---|---|
| 1 | `Plane.lean` | `Rigid` (rotation + translation as a structure; composition, inverse, action, isometry, `sda` kept), `rigidOfSegment` (the motion taking a segment onto a congruent one; uniqueness), triangles `Fin 3 → ℝ × ℝ` with `hull`, `interior`, open sides, corners; segment–triangle lemmas (the parameters of a segment inside a triangle form a closed interval; a segment entering the interior exits through one open side, transversally; a point on an open side is neither interior nor a corner; opposite sides of a line; the open wedge between two transversal lines at a point meets any neighbourhood); the chord bound `|PQ| ≥ |vQ| sin ∠PvQ`; angle additivity for collinear points (`cornerAngle` of `Triangle.lean`, by the law of cosines and Stewart's theorem) and the run bound (a straight segment through a fan of wedges at `v` subtends `< π`) | 1400 | planned |
| 2 | `Layout.lean` | the standard layout of every face of `s`, the unfolding motion `u_e` across each interior halfedge, cells `(σ, g)` with `tri`, `hull`, `across`; `ConsistentGerm`, `BoundaryGerm` and their transport by rigid motions; words, `devLength`, finiteness of words of bounded length (`Λ_L`), the constants `α_min`, `ℓ_min`, `N s L` | 600 | planned |
| 3 | `Overlay.lean` | `Overlay s t` (layouts of `t`'s faces with `t`'s lengths, cells, cover, `germ_int`, `germ_edge`, `germ_bdry`, `corner`, non-junk); the identity overlay of `s` on itself; **the walk**: along every halfedge of `t` a chain of cells, each entered through an open side and left through another, and the word it defines; `devLength word = t.ℓ x`; the fan bound `word.length ≤ N s (t.ℓ x)` | 1800 | planned |
| 4 | `OverlayFlip.lean` | `Overlay.flip`: a flip of a strictly convex diamond carries an overlay to an overlay of the flipped state (the two faces laid out as one convex quadrilateral `D` by `diamond_unique`; the new faces' layouts are the two triangles of the other diagonal; cells restricted to those overlapping the new face; every invariant re-proved by planar arguments: convexity of `D`, the diagonal is crossed by any segment from one half to the other, corners of `D` are never interior to a cell) | 3000 | planned |
| 5 | `Finite.lean` | `boundedLengthFinite : ∀ s, BoundedLengthFinite s`; `reach_finite`, unconditional `acc`, `flip_algorithm`, `exists_delaunay`, FlipQueue's `loop_acc`/`loop_result` without hypothesis | 300 | planned |
| 6 | docs | `DdgProofs.lean`, `Audit.lean`, README rows/table, CITED trims, survey status, formula tags | — | planned |

Goal 2 (Bobenko–Springborn uniqueness of the Delaunay tessellation / two Delaunay
triangulations differ by tie flips) is **not scheduled** within this layer's budget: on a surface
it needs either empty immersed disks (the surface) or the connectivity of the flip graph of all
geodesic triangulations with a given vertex set, neither of which this route provides. It stays
cited (`it-delaunay-existence-uniqueness`). If stages 1–5 finish early, the plan will be updated.

## Target statements (Lean)

```lean
-- Surface/Overlay.lean
structure Overlay (s t : MetricGluing F) : Type  -- see above
theorem Overlay.exists_word (O : Overlay s t) (x : Halfedge F) :
    ∃ w : Word F, w.numCross ≤ N s (t.ℓ x) ∧ devLength s w = t.ℓ x
def Overlay.id (s : MetricGluing F) : Overlay s s

-- Surface/OverlayFlip.lean
theorem Overlay.flip (O : Overlay s t) (h : ConvexFlip t t') : Nonempty (Overlay s t')

-- Surface/Finite.lean
theorem MetricGluing.boundedLengthFinite (s : MetricGluing F) : BoundedLengthFinite s
theorem MetricGluing.reach_finite (s : MetricGluing F) : (reach s).Finite
theorem MetricGluing.acc (s : MetricGluing F) : Acc (Function.swap DelaunayFlip) s
theorem MetricGluing.exists_delaunay (s : MetricGluing F) :
    ∃ t, ReflTransGen DelaunayFlip s t ∧ t.G.IsDelaunay t.ℓ
theorem MetricGluing.flip_algorithm' (s : MetricGluing F) : ... -- flip_algorithm without hfin
theorem loop_acc' / loop_result' -- FlipQueue without hfin
```

Hypotheses throughout: `[Fintype F] [DecidableEq F]`, `s : MetricGluing F` (a metric, so every
face is a nondegenerate triangle). Nothing else.

## What stays cited after this layer

* `it-delaunay-existence-uniqueness`: uniqueness of the Delaunay tessellation; the surface case
  of "Delaunay triangulations differ only by tie edges". Existence (surface, original vertex set,
  by flipping) becomes proved.
* `it-delaunay-flip-termination`: nothing (the surface case becomes proved), except the
  non-existence of a bound on the number of flips in `|V|`, which is a remark, not a theorem.

## Status log

* 2026-11-05: plan written; stage 1 started.
* 2026-11-05: stages 1–3 done and committed (`Plane`, `Segment`, `Angle`, `Quad`, `Count`, `Layout`,
  `Overlay`, `Fan`, `Relayout`: ~4,700 lines, no sorry). Stage 4 split into `FlipFace.lean` (face
  level, key crossing lemma) and `FlipEdge.lean` (edge level, assembly `Overlay.flip`); stage 5
  `Finite.lean` written against them. Deviation from the table: the planar file is four files,
  and the counting lemma got its own file.
