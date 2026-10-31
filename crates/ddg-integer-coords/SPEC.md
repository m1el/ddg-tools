# Integer Coordinates — contract

The contract between the Rust crate `crates/ddg-integer-coords`, the Lean oracle (`lean/oracle`),
the Lean proofs (`lean/proofs`), the Python binding (`crates/ddg-py`) and the Blender add-on for

> Mark Gillespie, Nicholas Sharp, Keenan Crane. *Integer Coordinates for Intrinsic Geometry
> Processing.* ACM Transactions on Graphics 40(6), Article 252, 2021.
> <https://doi.org/10.1145/3478513.3480522>
> (`../IntegerCoordinates/integer_coordinates.pdf`, "the paper", 13 pages plus the supplement:
> Appendix A in the paper, Appendices B–D, "Supp.", in the same PDF; equation and algorithm
> numbers are the paper's).

**Reference.** geometry-central at commit `ab116f1` (the commit the authors' demo pins; MIT,
Copyright (c) 2017-2019 Nicholas Sharp and the geometry-central contributors) in
`../IntegerCoordinates/intrinsic-triangulations-demo/deps/geometry-central/`, cited as
`gc: <file>:<line>`, and the demo itself (`../IntegerCoordinates/intrinsic-triangulations-demo`,
MIT, Copyright (c) 2021 Nicholas Sharp and Mark Gillespie), cited as `demo: <file>`. Both may be
read and ported; ported code keeps its notice in `LICENSES/`. The map of the reference is
`runs_tmp/ic-reference-map.md` (§11 summarizes where it departs from the paper).

This crate is the second representation of intrinsic triangulations in this workspace. The first,
signposts (Sharp, Soliman, Crane 2019), is `crates/ddg-intrinsic` with its contract
`crates/ddg-intrinsic/SPEC.md` ("IT-SPEC", decisions `N<n>`). Everything IT-SPEC says about the
input mesh, the halfedge conventions, Delaunay flipping, refinement with boundary and the overlay
holds here unless this file says otherwise. Decisions of this file have ids `I<n>`.

Conventions of `CLAUDE.md` apply: `f64` lengths, stored indices `u32`, explicit option structs, a
report from every operation that does a variable amount of work, bad input is an error naming the
element and never a panic, nothing repaired or normalized silently, deterministic results at any
thread count, formula tags `[formula: ic-…]` (§9).

## 1. The problem

`T⁰` is the input triangulation (a `ddg_core::TriMesh` with edge lengths; positions are never
needed). `T¹` is an intrinsic triangulation of the same polyhedral surface: a Δ-complex whose
vertices `V¹ ⊇ V⁰` are points of the surface, whose edges are geodesics and whose faces are flat
(paper Sec. 2.1). `V* = V¹ \ V⁰` are the inserted vertices.

The signpost representation stores the correspondence between `T⁰` and `T¹` as floating-point
directions, and recovers the connectivity of the common subdivision by tracing. Rounding can make a
trace miss its target (IT-SPEC N8, N25; the paper's Table 1: 96.0 % / 69.1 % success on Thingi10k
for signposts). The integer representation stores the correspondence as integers, from which the
**connectivity** of the common subdivision is recovered exactly; only the positions of the
crossings are computed in floating point (§3.2).

What must be true of this crate, in one sentence: every operation keeps the integer data equal to
the normal coordinates and roundabouts of *some* geodesic realization of `T⁰` over `T¹`, that
realization is the one the operation asked for wherever the operation is exact (flip, edge split
along a curve, removal), and within rounding of it where a floating-point decision is made (face
split, §4.4); the common subdivision extracted from the integers is a valid mesh with exactly the
connectivity those integers encode.

## 2. Data model

### 2.1 Connectivity

`T¹` is a halfedge Δ-complex with the conventions of IT-SPEC §3 and N9: halfedge `h = 3f + k` runs
from corner `k` of face `f` to corner `k + 1`, twins are stored, a triangle may use a vertex or an
edge twice, edges from a vertex to itself and vertices of degree one are allowed, oriented
manifold with or without boundary. **I1:** this connectivity and its topological primitives (flip,
face split, edge split, removal of a degree-3 vertex, with their index bookkeeping) are factored
out of `ddg-intrinsic` into one type that both representations embed (work package WP0); the
signpost crate's results must stay bitwise identical (its whole test suite and the fingerprints of
IT-SPEC N15).

**The shared interface** (a trait in `ddg-intrinsic`, implemented by both representations; what
`ddg-flip-geodesics`, `ddg-intrinsic-apps`, `ddg-uv` and the bindings program against; FlipOut
contract R1): connectivity, edge lengths, corner angles; `can_flip` / `flip_edge` with the
refusals of IT-SPEC N5; `split_edge` returning both halves; a constrained-edge mask (edges that
flips and Delaunay flipping leave alone, like boundary edges); `trace_edge_on_input -> EdgeTrace`
(the path of an intrinsic edge over `T⁰`); `vertex_location`. Guarantee: a flip keeps the stored
direction of every edge it does not flip.

**As built (WP0, `ddg-intrinsic` N50).** The connectivity is `ddg_intrinsic::Complex` (faces, twins,
halfedge↔edge maps, vertex→halfedge, the constrained mask); the interface is the trait
`ddg_intrinsic::Triangulation`. Every topological operation of `Complex` returns its index
bookkeeping: `flip_edge -> ComplexFlip` (refusal reasons from `flip_refusal(e)`, in the order
boundary, constrained, degree-one endpoint), `split_face -> FaceSplit`, `split_edge -> EdgeSplit`
(`halves = [e, e2]` from the old tail, both directed like the old edge), `remove_vertex ->
VertexRemoval` (interior, degree 3, no self-edge, no constrained spoke). **Removal renumbers:** the
merged face takes the smallest of the three face ids, and the last vertex, edges and faces move
into the freed ids; `VertexRemoval` lists every move (`apply_{halfedge,edge,face,vertex}_data`),
and this crate moves `n`, `r` and `q⁰` with it. The constrained mask lives in `Complex`; splitting
a constrained edge gives two constrained halves; flips refuse it with `FlipRefusal::Constrained`.
`Insertion::halves` carries the two halves of an edge split. The refusals and new lengths of flips
and splits are decided from `(&Complex, &[f64])` alone by `plan_flip`, `plan_face_split` and
`plan_edge_split`, which this crate calls, so both representations refuse the same operations and
compute bitwise the same lengths.

### 2.2 Stored quantities (paper Sec. 2)

| quantity | where | type | paper |
|---|---|---|---|
| `ℓ_ij > 0` | per edge of `T¹` | `f64` | Sec. 2.2 |
| `n_ij ≥ −1` | per edge of `T¹` | `i64` (I2) | Sec. 2.3 |
| `r_{aj} ∈ [0, deg₀(a))` | per halfedge of `T¹` whose tail `a ∈ V⁰` | `u32`, `NONE` for tails in `V*` | Sec. 2.4 |
| `q⁰_i` | per vertex (`V⁰`: `Vertex(v)`, so `vertex_location` is one lookup; `V*`: an `Edge` point with `t ∈ (0,1)` or a `Face` point with all barycentric coordinates in `(0,1)`) | `ddg_intrinsic::SurfacePoint` on `T⁰` | Sec. 2.2 |
| `T⁰` | shared, immutable | `Arc<ddg_intrinsic::InputSurface>` | — |

* `n_ij` is the number of transversal crossings of edge `ij ∈ T¹` with edges of `T⁰`, or `−1` if
  `ij` runs along an edge of `T⁰` (then no edge of `T⁰` crosses it). `n⁺ = max(n, 0)`,
  `n⁻ = −min(n, 0)` (`ic-normal-sign`). Note the paper's convention differs from Gillespie et al.
  2021 (`n = 0` for parallel edges); we follow this paper.
* **Local index** (Supp. B): the halfedges of `T⁰` around `a ∈ V⁰` are numbered counter-clockwise
  `0 … deg₀(a) − 1`. **I3:** the start is the vertex's `outgoing` halfedge in `InputSurface`; at a
  boundary vertex that is the boundary halfedge with the surface on its left (as IT-SPEC N2), the
  last index `deg₀ − 1` is the arriving boundary halfedge (as the reference), and the numbering does
  not wrap (the modulus is `0` at a boundary vertex, the oracle's `DEG = 0`).
* `r_{aj}` is the local index of the first halfedge of `T⁰` at or after `aj` counter-clockwise (it
  may be `aj` itself, when `aj` runs along an edge of `T⁰`).
* Initially `T¹ = T⁰`: same connectivity and lengths, `n ≡ −1`, `r_{ab}` = the local index of
  `ab`, `V* = ∅`.

### 2.3 Derived quantities

Corner coordinates and emanating counts at corner `k` of face `ijk` (Eq. 1, 2; Fig. 6):

```
c^{ij}_k = ½ ( max(0, n⁺_jk + n⁺_ki − n⁺_ij) − e^{jk}_i − e^{ki}_j )      [formula: ic-corner-coordinates]
e^{ij}_k = max(0, n⁺_ij − n⁺_jk − n⁺_ki)                                  [formula: ic-emanating]
```

`c^{ij}_k` counts arcs that cut corner `k`; `e^{ij}_k` counts arcs from vertex `k` to edge `ij`.
Combinatorial crossing `ζ = (ij→, p)`, `0 ≤ p < n_ij`, the `p`-th crossing along `ij` from `i`; its
reversal `ζ̄ = (ji→, n_ij − p − 1)` (`ic-reverse-crossing`). Geometric crossing
`z = (ij→, p, ab→, u, v)` with the barycentric coordinates `u` along `ab ∈ T⁰` and `v` along
`ij ∈ T¹` (Sec. 2.5). **Conventions (WP2):** `u` is measured in the input edge's own direction (as
`SurfacePoint::Edge`; the paper measures along the curve's direction, i.e. `u` or `1 − u`), `v`
along the crossed halfedge from its tail; `GeometricCrossing` also carries `sin_angle`, the sine of
the crossing angle, which says how well `u` and `v` are conditioned.

### 2.4 Validity (the invariant every operation keeps)

A state is **valid** when, in addition to IT-SPEC §5 items 1–2 (strict triangle inequalities,
angle sums):

1. every corner coordinate is an integer: the bracket of Eq. 1 is even (`ic-face-validity`).
   It is never negative, and is even exactly when `Σ n⁺` over the face is even or some corner
   emanates (Lean `faceValid_iff`);
2. in every face at most one corner has `e > 0` (two cevians from different corners cross). This
   follows from Eq. 2 (Lean `emanating_exclusive`) and is not checked separately;
3. an edge with `n = −1` has no crossings and both its halfedges with tails in `V⁰` have
   roundabouts pointing at the edge of `T⁰` it runs along;
4. every edge `ab ∈ T⁰` is extracted (Alg. 4) as a curve that starts at `a`, ends at `b`, passes
   only through vertices of `V*` that lie on it, and the roundabout labels agree at both ends;
   every combinatorial crossing belongs to exactly one such curve. In detail (as `certify` checks
   it, WP1):
   * around each `a ∈ V⁰`, counter-clockwise, the next roundabout is `r + n⁻ + e` of the current
     halfedge (`e` the emanating count of the corner at `a`); at an interior vertex mod `deg₀` with
     `Σ (n⁻ + e) = deg₀`; at a boundary vertex `r = 0` on the boundary halfedge and the chain ends
     at `deg₀ − 1`;
   * each input edge, followed from its tail, arrives at its head by exactly the end that the
     head's roundabouts give for the reverse direction;
   * an inserted vertex on an input edge has exactly two curve ends there and the vertices passed
     appear in the order of their parameter `t`; an inserted vertex inside a face has none;
   * every crossing and every `n = −1` edge is used exactly once.

Before the items, **ranges**: `n ≥ −1`, `n = −1` on boundary edges, no self-edge with `n = −1`; a
roundabout lies in `[0, deg₀)` exactly when its tail is in `V⁰`; `q⁰` of an inserted vertex as in
§2.2 (never a `Vertex` point, never a `Face` point with a zero coordinate); input vertices keep
their boundary status.

Items 1–3 are local and cheap; item 4 is the **certificate** (§8 layer 1). `certify()` checks all
four, in the order ranges, 1, 3, 4, and stops at the first failure. Its cost is linear,
about 100 ns per `|E| + Σ n⁺` (I4, measured in WP1: spot 2.1 ms, camelbad 9.8 ms; one certificate
≈ 1.7–1.9 |E| signpost flips). **Default (I4 decided):** certify at the end of every public
operation that does a variable amount of work (Delaunay flipping, refinement, before the common
subdivision, a FlipOut run), as a `*_certified` variant where the operation is the shared one
(WP3: `flip_to_delaunay_certified`; an `IntegerOptions::certify` switch does not exist yet); never per flip or
per split, also not in debug builds (it would make debug Delaunay flipping quadratic); tests
certify per operation on small meshes and per batch on large ones.

## 3. Reading the integers

### 3.1 Tracing (Supp. Alg. 1–4)

* `TraceFrom(ζ)` (Alg. 1, Fig. 7/20): enter the face across `ij→` at crossing `p`; with `k` the
  opposite vertex, the curve leaves through `ik` at `p` if `p < c^{kj}_i` … exactly the three cases
  of Alg. 1 (`ic-trace-step`). Terminates at a vertex.
* `ExtractCurve(ζ)` (Alg. 2): both directions, the combinatorial crossings, `IdentifyTracedEdge`,
  then the geometry by unfolding the strip and intersecting the straight segment with every
  interior edge (`ic-strip-intersection`). The strip unfolding reuses `ddg_intrinsic::formulas`
  (`it-unfold-diamond`, `it-ray-edge-intersection`).
* `IdentifyTracedEdge` (Alg. 3): walk back to a vertex of `V⁰`, then the local index
  `k = r_{ai} + (p − c^{ja}_i) + n⁻_{ai}` (`ic-identify-edge`), modulo `deg₀(a)` at an interior
  vertex. The halfedge is `Next(Next(ij)) = a→i` (Alg. 3's `aj` is a misnomer; the reference
  reads `roundabouts[he.next().next()]`, `ic:1983-1984`).
* `ExtractEdge(ab)` (Alg. 4): from the roundabouts at `a`, the halfedge of `T¹` just before `ab`,
  the crossing index `p = k + n_{ai} − r_{ai}`, then `ExtractCurve`, continued through vertices of
  `V*` that lie on `ab` (only edge splits create those, §4.5). The formula is correct as printed:
  it equals `c^{ja}_i + (k − r_{ai} − n⁻_{ai})`, since `c^{ja}_i = n⁺_{ai}` whenever `a` emanates;
  when the window wraps at an interior vertex use `k + deg₀`.

* `ExtractIntrinsicEdge(ij)` (the other direction, needed by FlipOut and by `trace_edge_on_input`;
  the reference has no such routine and builds the whole common subdivision instead, `ic:55-88`):
  for each of the `n⁺_ij` crossings of `ij`, the input edge crossing it is identified exactly from the
  integers (`ExtractCurve` through that crossing, Alg. 1 + 3), consecutive crossings are checked to
  share an input face, and the positions come from the unfolded input strip. The result is a
  polyline on `T⁰` whose connectivity is exact (`ic-extract-intrinsic-edge`). For `n_ij = −1` it is
  the input edge itself. **As built (WP2):** the walk from crossing `p` stops at a crossing of `ij`
  already identified, and crossings of `ij` it passes are kept in a sparse map, so the whole edge
  costs `O(n⁺_ij + walk steps)`; the first input face comes from the roundabout (the wedge between
  local indices `r − 1` and `r`) or from `q⁰`, the arrival is checked against the head's roundabout
  or `q⁰`; the length mismatch of the straight segment and the clamps are reported, never a
  failure; `EdgeTraceOptions::tolerance` is validated but not used, `max_crossings` is not used.

Every step works on integers and indices; a value outside its range (a crossing index `≥ n`, a
negative corner coordinate, a walk that revisits its start crossing without reaching a vertex) is
an error naming the face and the crossing, never a panic (`IntegerError::InvalidCoordinates`).
Per-edge operations must not compute the global bound `Σ n⁺ + |E¹|` (that makes all-edges work
quadratic): on valid faces Alg. 1's step is a bijection, so a walk either reaches a vertex or
returns to its start, which is detected.

### 3.2 What floating point decides

Only: the edge lengths after an operation (as in IT-SPEC), the positions of geometric crossings
(the parameters `u`, `v` of the strip intersection), the region of a face split (§4.4), the
circumcenter's location in refinement (§5), and the transfer matrices (§6). **I5:** a strip
intersection that falls outside `[0, 1]` by rounding is clamped and counted in the report
(`clamped_crossings`), never used to change connectivity (this covers the parameter along the
segment too; a halfedge parallel to the segment is cut at its middle and counted as a clamp of
`0.5`); with exact lengths the straight segment
lies inside the strip (the strip is the face sequence of the geodesic).

## 4. Operations

Every operation returns a report; a refused operation changes nothing (tested by comparing the
whole state).

### 4.1 Init

From a `TriMesh` + `Geometry` (as `IntrinsicTriangulation::new`): `T¹ = T⁰` (§2.2). Optional
intrinsic mollification first (§7, `ic-mollify`), reported.

### 4.2 Edge flip (Sec. 3.2, Eq. 3, 4)

Flippability as IT-SPEC §4 and N5 (strictly convex unfolded quadrilateral, both endpoints keep
degree ≥ 1, not a boundary edge). Then, for `ij` shared by `ijk`, `jil`, the new edge `kl`:

```
n_kl = c^{ji}_l + c^{ij}_k + ½|c^{il}_j − c^{ki}_j| + ½|c^{lj}_i − c^{jk}_i| − ½e^{ji}_l − ½e^{ij}_k
       + e^{lj}_i + e^{jk}_i + e^{il}_j + e^{ki}_j + n⁻_ij                     [formula: ic-flip-normal]
r_{kl} = (r_{ki} + e^{il}_k + n⁻_{ki}) mod deg₀(k)                           [formula: ic-roundabout-update]
```

(and `r_{lk}` symmetrically, only for tails in `V⁰`); `ℓ_kl` by `it-flip-length`. **Correction
to the paper:** its first two terms read `c^{jk}_l + c^{il}_k`, which name no corner of the two
faces; they are the corners at `l` of `jil` and at `k` of `ijk`, `c^{ji}_l + c^{ij}_k`, as the
reference computes (`nc:724-764`, `Clji`, `Ckij`) and Lean proves (`flip_count`). Eq. 3 can return
values below `−1` for combinatorially valid quadrilaterals with two or more `k`–`l` arcs (quad
`2,0,0,0,0` gives `−2`); straight arcs cannot do that, so the code treats `n_kl < −1` as
`InvalidCoordinates` naming the edge. The involution (`ic-flip-involution`) holds regardless. In code:
`plan_flip(&complex, &lengths, e)?` (refusal and `ℓ_kl`), the integers above from the pre-flip
state, then `complex.flip_edge(e)`, whose `ComplexFlip` says where `n_kl`, `r_kl`, `r_lk` go.
`can_flip_with(e, FlipShapeTest)` and `flip_edge_with(e, FlipRequest { shape, corner_angles })`
refuse by `plan_flip_with` (IT N52; FlipOut uses gc's relative test); every operation of this crate
keeps N5 and requests no angles. `Flip::new_corner_angles` is `Some` only when requested (from
`FlipPlan::new_corner_angles`, two `corner_angle`), so both representations return the same `Flip`
for the same request, and Delaunay flipping pays nothing for them (camelbad 9.47 → 8.89 ms). Exact: the
integers after the flip are those of the same realization over the new `T¹` (Lean tier C, §9).
Flipping back restores `n`, `r`, `ℓ` (integers exactly, lengths to rounding).

### 4.3 Delaunay flipping (Sec. 4.1)

IT-SPEC §4 and N6 unchanged (queue, tolerance, ties, limit), over this representation. On the same
input and the same queue order, the sequence of flips and the resulting `T¹` (connectivity and
lengths, bitwise: same formula) are identical whichever representation runs it. Tested against
`ddg-intrinsic`. Constrained edges (§2.1) are never flipped and are treated as Delaunay, as boundary
edges are. This is the provided method `Triangulation::flip_to_delaunay`, written once over the
trait (it is also `IntrinsicTriangulation::flip_to_delaunay`), so the identical flip sequence holds
by construction; `DelaunayReport::n_constrained` counts the constrained edges, which the other
statistics leave out. Only `FlipRefused` counts as a refused flip; any other error of a flip (an
`Overflow`, `InvalidCoordinates`) ends the flipping and is returned as
`IntrinsicError::Representation`, the failing flip changing nothing (55b0c3b). Measured (WP3,
load 3.3): Delaunay flipping 1.63 ms on spot (signposts 1.57), 9.29 ms on camelbad (8.87); random
flips 3.8 M/s and 3.2 M/s, as the signposts.

### 4.4 Face split (Sec. 3.3, Supp. Alg. 5, 6, App. A)

Insert `m` at barycentric `v` in `ijk ∈ F¹`. Gather the geometric crossings of the face's sides
(`ExtractGeometricCrossing`), find the region `R` of the face cut by the arcs that contains `m` by
line-side tests, and from it the new `n_{mi}, n_{mj}, n_{mk}` (Alg. 5 for `e = 0`, Alg. 6 for the
corner with `e > 0`; `ic-face-split-normal`), the roundabouts of new halfedges with tails in `V⁰`
(Eq. 4), the lengths from `plan_face_split` (`it-point-distances`, **I25**; the paper's Eq. 5,
`ic-tangent-length`, `‖w‖² = −ℓ²_ij w_i w_j − ℓ²_jk w_j w_k − ℓ²_ki w_k w_i`, is the same quantity
and is proved equal in Lean, but is not evaluated), and `q⁰_m` by `RecoverBarycentric` (App. A: minimum-norm `ξ` of the 3×ρ system of
Eq. 15, then Eq. 16; `ic-recover-barycentric`).

**Floating point.** The region test may pick a region next to the true one when `m` is within
rounding of an arc (Supp. Alg. 5 lines 10–16: keep the largest slack, zero the others). The result
is then valid (§2.4) and is the exact integer data of `m` placed in the chosen region; that
guarantee is proved (Lean tier B, §9), the choice itself is not. **I6:** the report records, per
split, the distance of `m` to the nearest arc in units of the face's shortest edge, and the split
is refused with `TooCloseToCurve` when it is below an option (default 0: never refused; set by the
callers that need the inserted point exactly, see §5).

**As built (WP4).** `split_face(f, bary, &FaceSplitOptions { too_close })` (and the trait's
`insert_vertex`; `locate_in_face` without inserting). Regions are enumerated from the integers
alone (`FaceArrangement`); the side tests are barycentric determinants in each chord's corner frame
(`ab − b·v₁ − a·v₂` for an arc around a corner, `β·v₁ − (1−β)·v₂` for a fan arc; computing them in
the plane picked a wrong region when an arc's ends were within `1e-12` of a corner). The region
with the largest margin wins, ties in the oracle's order, and the counts are `split_counts` of that
region (I15), so the result is valid whichever region is picked. Crossing parameters are clamped
to `[0, 1]` and counted (as I5). The I6 measure `curve_distance` is the distance of `m` to the
nearest bounding chord, in shortest sides. The region's input face is read from a gap of the region
on a side, or from the roundabout when all three sides run along input edges. Errors:
`TooCloseToCurve { face, distance, limit }`, `RecoverFailed { face, input_face, detail }`.
Measured: of 26 099 insertions (17 shapes, spot, camelbad, 7 scenarios), the integers equal
`from_intrinsic` of the signposts in all but 9, all within `1.2e-10` shortest sides of an arc
where rounding puts the point on different sides in the two representations (both valid);
the signposts snap places within `1e-12` of an input edge onto it (816 cases), so those states
differ by design; refusals identical; `q⁰` within `5.6e-10` mean edge lengths. Insertion costs
1–3.5 µs (signposts 0.9–1.7 µs; the integers run `extract_intrinsic_edge` on the three sides).

### 4.5 Edge split (Sec. 3.4, Alg. 7, Eq. 6)

`n_ij ≥ 0`: face split then flip (Fig. 10 top). `n_ij = −1`: explicit split along the shared
curve, `n_mj = n_mi = −1`, `n_mk = max(n_ki, n_jk, 0)` (`ic-split-edge-normal`; the `0` matters when both are `−1`), and the same on the
other side when interior; `q⁰_m = u_i q⁰_i + u_j q⁰_j` along the input edge. Boundary edges are
always of the second kind. Lengths and refusal from `plan_edge_split`, connectivity from
`Complex::split_edge` (I20: also for `n_ij ≥ 0`, without the face-split-then-flip route).
**`n_ij ≥ 0`, as built (WP3):** the split needs only `before`, the number of crossings of `ij`
before the split point (those of `extract_intrinsic_edge(ij)` with `distance/ℓ < t`), and `q⁰_m`,
interpolated between the two neighbouring crossings inside their input face. The halves get
`before` and `n_ij − before`; the new edge to the opposite corner `c` of each side
(`ic-split-crossed-edge`; Eq. 6 is its `n = −1` case, tested):

```
n_mc = c_c + e_a + e_b + max(0, c_a − before) + max(0, before − (n_ab − c_b))      [formula: ic-split-crossed-edge]
```

with `a b` the split edge, `c_*`, `e_*` the corner and emanating counts of that side's face. A
split point within rounding of a crossing is refused (`InsertionRefused`), nothing changed.

### 4.6 Vertex removal (Sec. 3.5, Alg. 8, Thm. D.1)

Only `i ∈ V*`. Flip flippable edges at `i` until its degree is 3 (interior) or the boundary
analogue (I7: the reference's rule for boundary vertices, `runs_tmp/ic-reference-map.md`), then
delete. No integer data changes beyond the flips. Thm. D.1 proves that a simplicial star with
angle sum `2π` and degree `≥ 4` has a spoke whose quadrilateral has all angles `≤ π`; it need not
be strictly convex (Lean `cross_star_no_strict`: in the degree-4 cross every inner corner is `π`),
so N5's strict flippability can fail. At an interior vertex made by an edge split (a straight
star) the two spokes at angle `π` are never flippable; the others may be (**I27**).
For a star that is not simplicial the paper conjectures it (open, `ic-vertex-removal-delta`).
**I8:** the flip choice is deterministic (I24: the reference's score, self-loops first, ties by
edge index; only flips that lower the degree count; the reference's `checkFlip` gives refused,
fixed or degree-one edges the score `+∞`), the loop is bounded by `10 deg(i)` flips, and a star
where no such flip exists is refused with `RemovalStuck { vertex, degree }`, the state restored
(by a clone of `Complex`, `O(|E|)`; a journal of inverse flips does not restore `vertex_halfedge`).
Counted over every removal of the evaluation (§10). WP3: 3884 removals with 2.4 flips each (at most
14); of face-inserted vertices 1148 removed, 12 stuck.

### 4.7 Moving an inserted vertex (Sec. 3.6)

`x = exp_i(v)` traced on `T¹` (IT-SPEC tracing over the intrinsic triangulation needs signposts;
here **I9:** the tangent vector is laid out in the faces around `i` from the edge lengths, i.e.
the trace is `ddg_intrinsic::AngleFrames::of(&t).trace(start, dir, len, &opts)`, the same kernel
with signposts computed from corner angles. Frame convention: direction 0 at a vertex is
`Complex::outgoing_halfedge`, the others the running sum of corner angles times `2π/Θ` (`π/Θ` at
the boundary), `Θ` from the lengths; a frame is valid only for the current connectivity, since a
flip at an interior vertex can move its outgoing halfedge. Measured against stored signposts:
within 6.2e-12 rad on 17 shapes, traces end within 6.2e-12 of the mean edge length), then `SplitFace(x)`, then `RemoveVertex(i)`. The optimal Delaunay
iteration of IT-SPEC N21 then runs on this representation unchanged (not ported yet).

**As built (WP4).** `move_vertex(v, displacement, &FaceSplitOptions) -> VertexMoveReport`: the
displacement is given in the `AngleFrames` frame of `v` (direction wrapped into `[0, 2π)`); a path
that ends on an edge is an edge split; the moved vertex keeps its id (removal renumbers the last
vertex into the freed one); boundary vertices are `Unsupported` (I7); a refused move
(`RemovalStuck`, a path leaving the surface) restores the whole state from a copy. Measured: 1603
moves inside the star agree with the signposts' `move_vertex` within `7.6e-13` mean edge lengths;
951 moves beyond the star (which the signposts cannot do) certify.

### 4.8 Common subdivision (Sec. 3.7, Supp. Alg. 9)

Per face of `T¹`, oriented so that `n_ij ≥ n_jk, n_ki` (Alg. 9 line 5): case 1 (`e^{ij}_k = 0`,
corner-to-corner pairing) or case 2 (fan from the one corner with `e > 0`), Fig. 11. Vertex
positions: `V¹` from `q⁰` and the vertex itself, crossings from `ExtractEdge` over every `ab ∈ E⁰`
(Alg. 9 lines 15–20). Output: `ddg_intrinsic::Overlay` (the same type, so IT-SPEC N12's tests and
every consumer apply unchanged), plus the interpolation matrices `P₀` (`|V^S| × |V⁰|`) and `P₁`
(`|V^S| × |V¹|`), rows of 1 to 3 entries summing to one (Sec. 4.3). A face whose arcs fail the
pairing is impossible for a valid state; if it happens it is an error naming the face (it would
mean §2.4 was violated, which `certify` must have caught).

**As built (WP5).** `common_subdivision() -> CommonSubdivision { overlay, fingerprint, report }`;
`interpolation_matrices(&CommonSubdivision) -> (P₀, P₁)` as `ddg_core::SparseMatrix`, refused when
the fingerprint is not the current state's. No rotation to make `n_ij` largest is needed: the
pairing depends only on `c` and `e` (`face_arcs`, `ic-face-decoding`). The input face of a piece
comes from the orientation of its crossings (the tail side of a `T¹` halfedge lies in the
right-hand input face, the head side in the left-hand one), with a fallback through
`extract_intrinsic_edge` for faces whose sides all have `n⁺ = 0`. Overlay vertices: `V¹` first
(so the first `|V¹|` rows of `P₁` are the identity), then the crossings edge by edge from the tail
of the stored halfedge, which are the signpost overlay's ids. A self-edge crossing collapses to one
entry of `P₁`. Measured on the 266 states: connectivity identical to the signposts' overlay as
labelled polygons, positions `|Δ| · sin θ ≤ 3.5e-10`, the N12 checks pass; 24 ms on camelbad
against 8.7 ms for the signposts' traces plus overlay (the integers also pay for Alg. 2's strip
geometry).

### 4.9 Point queries

`T¹ → T⁰`: a point of a face of `T¹` → the region of the face → `RecoverBarycentric` (as §4.4
without inserting). `T⁰ → T¹`: by the common subdivision of the containing input face (I10:
lazily, per input face). `ExtractEdge` of its three edges alone cannot do it, since `T¹` edges
between inserted vertices inside the face are invisible from them; as built (WP5) a search over
`T¹` faces starts at the face holding the start of `ExtractEdge(3g)`, uses
`extract_intrinsic_edge` of each face's sides, and crosses every side of a piece that is not part
of an input edge (`PointLocator` caches the pieces per input face). Results are snapped:
a barycentric coordinate `≤ LOCATE_SNAP = 1e-12` makes an edge or vertex point. Same signatures as
`ddg-intrinsic`'s `point_on_input` / `point_on_triangulation`.

## 5. Delaunay refinement (Sec. 4.2, Supp. Alg. 10, App. C)

Chew's second algorithm, intrinsic, with boundary:

1. Flip to Delaunay (§4.3).
2. Pick a triangle below the angle bound `θ_min` (I11: the order of IT-SPEC's refinement queue, so
   that both representations insert the same points on the same input where the region tests agree).
3. Circumcenter in homogeneous barycentric coordinates `v̂_i = ℓ²_jk(ℓ²_ij + ℓ²_ki − ℓ²_jk)` and
   normalized (Eq. 9, 10; `it-circumcenter`); the tangent vector from the barycenter; `exp` from the
   barycenter (Fig. 13 right; tracing as I9).
4. Inside the surface: insert it (§4.4). Hits a boundary edge `lm` first: split `lm` at its midpoint
   (§4.5), flip to Delaunay, then remove every inserted interior vertex within Dijkstra distance
   `ℓ_lm` (edge-graph distance; Xia 2013 Cor. 1: the geodesic ball of radius `ℓ_lm/2` lies inside
   the Dijkstra ball of radius `ℓ_lm`) of the new vertex (§4.6), flip to Delaunay.
5. Narrow cones: triangles incident to exactly one vertex of angle sum below 60° or contained in an
   input triangle incident to one are not refined and are excluded from the reported minimum angle
   (Sec. 4.2 last paragraph); IT-SPEC N4 records the signpost crate's rule; I12 says which rule
   this crate uses after comparing both with the reference.

Constrained edges are treated like boundary edges: encroachment splits them, circumcenter traces
stop at them (later, for refinement that keeps FlipOut paths as edges; FlipOut contract F9).

Options and report as IT-SPEC (`RefineOptions`, `RefineReport`): angle bound (paper 25° in the
evaluation; Thm. C.1 gives 30° for closed surfaces with cone angles ≥ 60°), optional circumradius
bound, insertion limit (the reference's default `−10` = `10 |V⁰|`), and additionally the removals,
refused removals (I8), clamped crossings (I5) and the certificate's result.

**As built (WP6, 49bac6c, 91ca49b).** One driver for both representations (I11):
`ddg_intrinsic::refine::refine<T: Triangulation>` with `RefineRules { Ruppert, Chew }`.
- `Ruppert` is IT-SPEC N1/N4/N14 unchanged and stays the signposts' default (bitwise pins pass).
- `Chew` is this section and the integer default (`IntegerTriangulation::refine`,
  `IntegerRefineOptions`, certificate at the end). It splits at the midpoint when the circumcenter
  path leaves through a boundary edge, then Delaunay, then removes the inserted interior vertices at
  Dijkstra distance `< ℓ` (strict, as Alg. 10 line 20; the reference uses `≤`), then Delaunay.
- A circumcenter on an interior edge splits it at the foot without removal (the reference removes
  after every split).
- `RecoverFailed` and `TooCloseToCurve` are insertion refusals: the face is skipped, not an error.
- Trace failures are counted under both rules.
- Trait additions (provided methods): `trace_geodesic`, `remove_inserted_vertex`,
  `refine_with_rules`. The I24 removal choice is shared (`ddg_intrinsic::remove::removal_flip_choice`),
  so both representations remove by identical flips.

Measured over 256 differential runs (17 generator shapes, spot, camelbad, small-angle cones,
bob, blub, hand, …; both rules; 25° and 30°):
- 255 made the same operations with the same connectivity, bitwise-equal lengths and counts. Of
  these, 242 equal `from_intrinsic` of the signposts exactly; 13 differ only by the signposts'
  snaps onto input edges.
- The one divergence is a circumcenter exactly on an interior edge at its crossing with an input
  edge (I28).
- 25° is reached wherever possible, and 30° on closed shapes with all angle sums `≥ 60°` (Thm. C.1).
- 2307 removals, 175 refused (mostly I27 straight stars); 0 trace failures.
- Integer refinement costs 0.93–1.25× the signposts' (the paper: about 2.7×); the certificate adds
  about 20–80 %.

## 6. Attribute transfer (Sec. 4.3, Eq. 11–14)

Given `f` on `V¹` (piecewise linear on `T¹`), `f̂` on `V⁰` minimizing
`‖f − f̂‖²_{L²} = (P₁f − P₀f̂)ᵀ M_S (P₁f − P₀f̂)` (Eq. 13) solves `P₀ᵀ M_S P₀ f̂ = P₀ᵀ M_S P₁ f`
(Eq. 14), `M_S` the Galerkin mass matrix of the common subdivision (`ddg_intrinsic_apps::fem`
conventions). The paper writes `M₀ f̂ = P₀ᵀ M_S P₁ f`; `P₀ᵀ M_S P₀ = M₀` holds exactly because the
pieces of an input triangle are linear on it (`ic-transfer-mass-identity`: the finite identity
for the pieces of Alg. 9, proved; `ic-transfer-mass-integral`: the tiling and the P1 integral,
cited, Ciarlet 1978), and is checked numerically. Also the reverse direction (`T⁰ → T¹`) and plain copying, for the comparison of
Fig. 14. The factorization of `M₀` is kept for repeated transfers. Report: residual, the size of
`S`.

**As built (WP7, b463c96).** `Transfer::new(&tri, &CommonSubdivision)` (or `from_state`) builds
`M_S`, `P₀`, `P₁`, `M₀ = P₀ᵀ M_S P₀`, `M₁ = P₁ᵀ M_S P₁` and Cholesky-factorizes `M₀`, `M₁` once;
`to_input` (Eq. 14) and `to_intrinsic` (`M₁ ĝ = P₁ᵀ M_S P₀ g`), `copy_to_input`,
`copy_to_intrinsic`; a stale subdivision (fingerprint) is refused. `M_S` is unnormalized
(`Σ M_S` = total area; the area normalization of CLAUDE.md belongs to the cones problem); piece
areas are `|det β| · Area(T)` with `β` the barycentric rows of the fan triangle in its input face,
so no positions are needed. Measured on 133 states: `P₀ᵀ M_S P₀` equals the input mesh's own
Galerkin mass matrix to `2.1e-15` (relative to the larger diagonal entry), piece areas sum to the
total area within `1.1e-14`; constants reproduced to `5e-15` both ways, normal-equation residual
`≤ 2.6e-16`. Building everything on camelbad (`|S|` 62 825) takes about 250 ms, one transfer
2.8 ms.

## 7. Robustness preprocessing (Sec. 3.9)

Intrinsic mollification (Sharp and Crane 2020b) with `ε = 10⁻⁵ h̄` (`h̄` the mean edge length) by
default, as in the evaluation; uses `ddg_core::mollify::mollify_edge_lengths` (I13: moved there in
WP0, errors as `MollifyError`; `ddg_bff` keeps a wrapper). The lengths it changes are reported; the
operation is off by default in the library (nothing silently changed) and on in the evaluation and
in the Python/Blender defaults, as the paper uses it.

## 8. Verification layers

1. **Certificate** (§2.4): `certify()` after every operation in the tests and in debug builds; the
   evaluation runs it on every output.
2. **Differential test against signposts**: on well-conditioned meshes, the same flips give the
   same `T¹`; the common subdivision from integers and the overlay from `ddg-intrinsic`'s traces
   have the same connectivity (as labelled polygons) and positions within `1e-9` of the total
   area's square root; the point queries agree. Crossing positions are compared conditioned:
   `|Δ| · sin θ ≤ 1e-9` with `θ` the crossing angle (WP2: raw differences above `1e-9` occur only at
   grazing crossings, `sin θ ~ 1e-4`, where the signposts' `~1e-12` rad angle error is amplified;
   worst conditioned difference `1.8e-10` against signposts, `3.5e-10` between the two directions
   and against the overlay, over 272 651 crossings on 266 states).
3. **Brute-force recount**: for small inputs (`< 500` faces), after random sequences of flips,
   splits and removals, recompute every `n_ij` from scratch by tracing every input edge in exact
   rational arithmetic over the unfolded strips (only where the lengths are rational, i.e. on the
   generators' integer-coordinate grids) and compare with the stored integers.
4. **Lean oracle** (`lean/oracle`, no Mathlib): an integer implementation of Eq. 1–4, 6, Alg. 1–3
   and Alg. 9's per-face pairing, written from the paper; Rust is fuzzed against it with exact
   equality (no tolerance: integers).
5. **Lean proofs** (`lean/proofs`, Mathlib): §9.
6. **Against the reference** (§10.3).

## 9. Lean: the combinatorial model and what is proved

### 9.1 The model

A **face model**: a triangle with sides `0, 1, 2`, and a finite set of pairwise disjoint *arcs*,
each arc a pair of endpoints where an endpoint is a corner or a position on a side (sides carry a
linear order of their crossing points). Disjointness is combinatorial: two arcs are disjoint iff
their endpoints do not interleave on the boundary cycle of the triangle (the chord-interleaving
lemma below makes this the geometric statement for straight arcs). Normal coordinates of a face
model: the number of arc endpoints on each side.

A **complex model**: a finite Δ-complex given by gluing data (faces, the side pairing with
orientation reversal), one face model per face, with matching crossing counts on glued sides
(crossing `p` of a side is crossing `n − p − 1` of its twin). Curves: the maximal chains of arcs
glued through sides; their ends are corners. Labels: at each vertex of `V⁰`, the cyclic order of
`T⁰` halfedges and the roundabout map.

### 9.2 Theorems (ids, and the tier of the effort estimate)

| id | statement | tier |
|---|---|---|
| `ic-normal-sign` | `n⁺ − n⁻ = n`, `n⁺ n⁻ = 0` | A |
| `ic-corner-coordinates`, `ic-emanating` | for a face model, Eq. 1–2 compute exactly the arc counts by type; they are non-negative integers | A/B |
| `ic-face-decoding` | the three side counts (with the "at most one corner emanates" condition) determine the face model up to the order-preserving relabelling: existence and uniqueness (Fig. 11) | B |
| `ic-chord-interleave` | two straight chords of a convex polygon cross iff their endpoints interleave on its boundary (Mathlib, Euclidean plane) | B |
| `ic-cevians-cross` | segments from two different corners of a triangle to their opposite sides cross | B |
| `ic-reverse-crossing` | `ζ̄̄ = ζ`, and the gluing respects it | A |
| `ic-trace-step` | Alg. 1's three cases follow the unique arc of the face model through crossing `p` | B |
| `ic-trace-terminates` | in a complex model, `TraceFrom` ends at a corner after at most `Σ n⁺` steps and visits each crossing at most once | C (WP9; `CITED.md` lists it as pending until then) |
| `ic-identify-edge` | Alg. 3's index is the local index of the curve's `T⁰` halfedge | C (WP9; `CITED.md` lists it as pending until then; tested against signpost traces) |
| `ic-flip-normal` | in a convex quadrilateral (two glued face models, straight arcs), the arc system re-cut along the other diagonal has the counts of Eq. 3 (with the corrected first terms, §4.2) | A (identity) and the count in a quadrilateral model `Quad` (`flip_count`): proved; gluing two face models and re-cutting: C |
| `ic-flip-involution` | Eq. 3 applied twice returns `n_ij`, for counts coming from a face-model pair | A |
| `ic-roundabout-update` | Eq. 4 gives the local index after the flip | C |
| `ic-face-split-normal` | Alg. 5/6: inserting a vertex in any region `R` of a face model and joining it to the corners by arcs disjoint from the others gives the stated counts; hence the result is valid whichever region floating point picks | B |
| `ic-split-edge-normal` | Eq. 6 | A |
| `ic-split-crossed-edge` | §4.5: splitting a crossed edge after `before` crossings gives `n_mc` as stated, for any `0 ≤ before ≤ n_ab` of a valid face (Eq. 6 the `n = −1` case) | B (from the face model and chord interleaving; pending) |
| `ic-strip-intersection` | the straight segment through an unfolded strip meets each interior edge once, at the computed parameters (both directions, §3.1) | A (pending) |
| `ic-extract-intrinsic-edge` | `ExtractIntrinsicEdge` returns the input edges crossed by the intrinsic edge, in order | C (pending, WP9) |
| `ic-recover-barycentric` | Eq. 15–16 reproduce `q⁰` exactly for any solution `ξ`, when the corner positions are exact (affine combination) | A |
| `ic-tangent-length` | Eq. 5 (law of cosines in barycentric form) | A (proved equal to `it-point-distances`, which the code evaluates; I25) |
| `ic-transfer-mass-identity`, `ic-transfer-normal-equations` | the finite identity behind `P₀ᵀ M_S P₀ = M₀`; Eq. 14 is the minimizer of Eq. 13 | A |
| `ic-transfer-mass-integral` | the pieces tile the input triangle and the P1 mass integral | CITED (Ciarlet 1978) |
| `ic-vertex-removal` | Thm. D.1: a simplicial star of angle sum `2π` and degree `> 3` has a spoke whose quadrilateral has all angles `≤ π` (strictly convex when no inner corner is `π`; the strict version is false in general) | B′ |
| `ic-vertex-removal-delta` | the Δ-complex case (paper: conjecture) | attempt, bounded; else CITED as open with I8's runtime refusal |
| `ic-normal-isotopy` | complex models ↔ normal isotopy classes of curves on the surface (Kneser 1929, Haken 1961; Schaefer, Sedgwick, Štefankovič 2002; Erickson, Nayyeri 2013) | D, CITED |
| `ic-refinement-termination` | Thm. C.1 with Lemmas C.2–C.5; Chew 1993, Shewchuk 1997; boundary unproved in the paper | E, CITED |
| `ic-dijkstra-ball` | Xia 2013 Cor. 1 | E, CITED |
| `ic-mollify` | Sharp, Crane 2020b | CITED; the formula itself is proved as `intrinsic-mollification` |

### 9.3 Existing formalizations (survey of 2026-10-31, `runs_tmp/ic-formalization-survey.md`)

Nothing in §9.2 is formalized in any proof assistant found (Lean/Mathlib, Isabelle AFP, Coq,
HOL Light, Mizar), except as noted; the 2026 items of the survey came from search summaries and are
re-checked before they are cited.

* Mathlib (Apache-2.0) has the blocks, not the theorems: `Sbtw`, `SOppSide`, `Cospherical`,
  circumcentre, orthogonal projection and Gram matrices (tier A transfer items become cheap), a
  definitions-only `Geometry/Polygon/Basic.lean`; no Jordan curve theorem, no Δ-complexes or
  halfedge structures, no chord-crossing lemma. `ic-chord-interleave`, `ic-cevians-cross`,
  `ic-vertex-removal` are written from scratch by coordinates (EmptyHexagonLean's orientation
  theory read first for structure; its licence is unknown, so nothing is copied).
* Tier C's model: no Lean halfedge/hypermap library. The structural blueprint is Dufourd and
  Bertot's Coq hypermap development (ITP 2010: flips preserve the hypermap, triangulation and
  embedding invariants; planar Lawson termination via the lifted volume); also Gonthier's
  four-colour hypermaps (Coq) and `mccorvie/classification-of-surfaces` (Lean 4, Apache-2.0) for
  dart-based surface data. Structure may be ported, text is not copied.
* Tier D (normal curves ↔ isotopy), tier E (Chew/Ruppert/Shewchuk termination, Xia's stretch
  factor, Bobenko–Springborn existence/uniqueness, Rippa), FlipOut's guarantees: no formalization
  anywhere; they stay cited. Planar Lawson termination (Coq) is added to the citation of
  `it-delaunay-flip-termination` as a cross-reference, not as a proof of the surface statement.
* If a Jordan curve theorem is ever needed: `rkirov/jordan_pick` (Lean 4, Apache-2.0, polygonal and
  continuous JCT, claims sorry-free); tier C as designed (combinatorial disjointness by
  interleaving) does not need it.

Tiers: A formulas, B one face, B′ Thm. D.1, C the complex model, D and E cited. Work packages WP8
(A, B, B′) and WP9 (C) in §12; WP9 is last and can be cut without blocking anything else.

## 10. Measurements

### 10.1 What the paper reports (the floor to match or beat)

* Single core, Intel i9-9980XE: refinement of the Fig. 15 chair 0.2 s, of Fig. 16 (top) 0.6 s;
  Thingi10k 719791 (Fig. 19): 4 min with integer coordinates, 1.5 min with signposts.
* Table 1, Thingi10k, 7696 valid manifold meshes, mollified `10⁻⁵`: intrinsic Delaunay
  triangulation and refinement to 25°, success = the common subdivision's connectivity recovered:
  explicit overlay 100 % / –, signposts 96.0 % / 69.1 %, integer coordinates 100 % / 100 % (one
  refinement's subdivision, ~30 M vertices, not extractable for memory, Fig. 19). Signposts
  evaluated without boundary. **Reproduced (WP10):** integers 99.97 % / 99.97 % of our 7745
  meshes; the two misses are non-terminating refinements in the shared driver (I30), not
  representation failures (§10.2).
* Size: the common subdivision after refinement has mean `|V^S| / |V⁰|` 20, 95th percentile 45;
  the refined `T¹` has 3.7 and 7.8.
* Fig. 14: Poisson `Δf = g`, `g = sin(πx) sin(πy)` on the unit square, random low-quality meshes
  from random edge splits, `2⁶ … 2¹¹` vertices, 100 trials each: error of the solution on the
  original mesh vs on the IDT copied back vs the IDT transferred by Eq. 14; the IDT is 11.5× and
  the `L²` transfer 28.1× more accurate than the original at the marked size.
* Sec. 5.2: FlipOut geodesics on the integer representation, ~62 000 trials on Thingi10k, all
  connected. Measured in the FlipOut contract (`crates/ddg-flip-geodesics/SPEC.md` §9), §13 Q4.

### 10.2 Our plan

* **Thingi10k evaluation** (§13 Q1): both representations of this workspace, the reference, all
  manifold meshes, mollified, IDT and refinement to 25°: success (certificate passes and the common
  subdivision is a valid mesh, IT-SPEC N12's checks), time, sizes, removals refused (I8), clamped
  crossings (I5). Runs in parallel, one mesh per job, a 10-minute cap per job, logs per job; the
  timings of a 50-mesh subset alone, single thread.
  **Result (WP10, 307308a–a9fa390; tables in the crate README, runs in `runs_tmp/wp10/`).**
  7745 meshes, mollified `1e-5`, 25°, no insertion limit (the paper's benchmark setting; with 10
  per vertex, 147 meshes hit the limit). Success = common subdivision recovered and N12's checks:

  | | IDT | refinement |
  |---|---|---|
  | integers (default Delaunay test / the authors') | 99.97 % / 99.97 % | 99.97 % / 99.97 % |
  | signposts, Ruppert / Chew rules (with boundary) | 100.00 % / 99.97 % | 99.78 % / 99.78 % |
  | demo, integer backend (7444 meshes ≤ 200k faces) | 100.00 % | 99.97 % |
  | demo, signpost backend (same set, not mollified) | 87.86 % | 20.71 % |

  Integer failures: only 176193 and 190253, 10-minute time-outs of a refinement that does not
  terminate (I30). No certificate failure, invalid subdivision, trace failure, overflow or
  out-of-memory; 100 % on the 406 meshes with boundary. Signpost Ruppert failures: 15 trace
  failures plus the same two. Sizes after refinement: `|V^S|/|V⁰|` mean 19.83, p95 46.37 (paper 20,
  45); `|V¹|/|V⁰|` 3.80, 8.00 (paper 3.7, 7.8). Counters over the integer runs: 140 M insertions,
  184 706 removals, 1266 removals refused (110 meshes), 2840 clamped crossings in refinement and
  10 353 in the refined subdivisions (0 in the IDT ones), 791 moves (I16), 7 face splits refused
  (`RecoverFailed`, meshes 551019 and 97660; 551019 ends at 21.17° with the default Delaunay test
  and reaches 25° with the authors'). 169 meshes end with corners below 25° at narrow cones;
  Chew's rules leave 16 059 diametrally encroached boundary edges on 276 of the 406 boundary
  meshes.
  Single core, alone (`taskset`, outside load median 1.03): Fig. 15 0.029 s, Fig. 16 0.096 s,
  719791 58.6 s (paper 0.2 s, 0.6 s, 4 min); signposts 0.021, 0.064, 36.1 s. 719791's common
  subdivision (33.9 M vertices, the one the paper could not extract) extracts in 16.1 s at 12 GB
  peak. On the 50-mesh subset our integer refinement costs 1.7× the signposts' (paper 2.7×) and
  the demo's integer refinement a median 2.8× ours. Performance outlier: mesh 54621 (1450 boundary
  edges, 12 951 Chew removals) takes 28 s with integers and 65 s with Chew signposts against 0.6 s
  Ruppert and 2.1 s in the demo; the removal follows Alg. 10 line 20 (Dijkstra ball `< ℓ` after the
  insertion loop) where the reference removes after every split with `≤`.
* **Fig. 14** reproduced exactly as described (§13 has no question: the generator is ours, the
  random seeds are listed).
  **Result (WP7, 406c0e0, `examples/fig14.rs`).** Sizes `round(2^(k/2))`, `k = 12 … 22` (the
  figure's half-octave steps), 100 trials each, seeds `1000k + t`. Generator (the paper gives none,
  and the demo's `benchmark/` has only the Thingi10k run): the square cut by a diagonal, then
  repeated splits of the nearest edge of the face containing a uniform random point, at its
  projection clamped to `[0.05, 0.95]`. Absolute `L²` error against `f = −g/(2π²)`, by degree-4
  quadrature on the common subdivision; means over trials:

  | `|V|` | original | IDT copied back | IDT, `L²` transfer | IDT itself | orig / copy | orig / `L²` |
  |---|---|---|---|---|---|---|
  | 64 | 9.36e-3 | 5.18e-3 | 3.89e-3 | 3.13e-3 | 1.81 | 2.41 |
  | 181 (marked, `2^7.5`) | 6.93e-3 | 2.54e-3 | 1.70e-3 | 1.15e-3 | 2.73 | 4.09 |
  | 512 | 5.90e-3 | 1.34e-3 | 8.36e-4 | 4.13e-4 | 4.40 | 7.06 |
  | 2048 | 5.51e-3 | 6.30e-4 | 3.89e-4 | 1.02e-4 | 8.76 | 14.18 |

  The trends are the paper's (the original's error nearly flat, the IDT's falling like `1/|V|`,
  the transfer better than copying by 1.3–1.6 at every size, never better than the IDT itself);
  the factors at the marked size are 2.7× and 4.1× instead of 11.5× and 28.1×. The paper's
  original-mesh error is about `3e-2` (its meshes show long fans from the corners), ours about
  `7e-3`: our generator makes less degenerate meshes, and the ratios against the original depend
  on it. Not a defect of the transfer; matching the factors would need guessing the generator.
* **Benches**: flip, Delaunay flipping, refinement, common subdivision, transfer, on `data/models`,
  `data/intrinsic-reference` and the Fig. 15/16/19 Thingi models, single thread alone; the same
  benches for `ddg-intrinsic` side by side.

### 10.3 Against the reference

The demo with `--backend integer --noGUI` writes `faceInds.dmat`, `faceLengths.dmat`,
`vertexPositions.dmat`, `laplace.spmat`, `interpolate.spmat` (the README says `.mat`),
`{InputToIntrinsic,IntrinsicToInput}_{lhs,rhs}.spmat`, `common_subdivision.obj`, `stats.tsv`; dense
files start `# dense R C`, sparse ones `# sparse R C` with 1-indexed triplets (map §6). It cannot
dump `n` or `r`. Compare on the same inputs (both mollified with `1e-5`):

* IDT: edge lengths as sets (IT-SPEC §5 item 8) and, since `|V^S| = |V¹| + Σ n⁺`, the sum of the
  normal coordinates; the common subdivision's `V`, `E`, `F`.
* Refinement: vertex counts, minimum angle over the valid faces (`stats.tsv`), subdivision sizes
  (the orders differ, I23: statistics, not equality).
* Transfer: the demo assembles `M_S` from extrinsic positions, we from lengths, so the systems
  agree to the rounding of those two (`1e-9` relative), not bitwise.
* The reference's own tests (map §6) are ported as tests here, on its assets `fox.ply` (closed) and
  `cat_head.obj` (12 boundary edges), copied with the MIT notice.

Built outside this repository (polyscope and happly fetched at the demo's pinned commits,
2026-10-31).

**Result (WP10).** Two differences of the demo matter for any comparison: its signpost backend
does not mollify (its integer backend does), and it flips when the cotan weight is below `−1e-6`
(`triangleTestEPS`), so IDT parity needs `DelaunayCriterion::REFERENCE` (the `integer-ref`
backend of `thingi_eval`). With it, on 6004 meshes up to 20k faces: `|V¹|` equal on all; `|V^S|`
and the subdivision's `V`, `E`, `F` equal on 4899 (81.6 %), and on 1392 of the 1393 meshes without
near-ties; the rest differ by the direction of near-tie flips (I23). Edge-length sets match to
`1e-9` on 4581 meshes. Refinement (7440 meshes): `|V¹|` ours/demo mean 0.9998, median 1.0000,
within 1 % on 6531; `Σ|V^S|` −0.14 %; the demo ends below 25° on 9 meshes, ours on none. The demo's
own failures are its assertions ("could not identify source face" 5336, "must sum to an even
number" 1213), 40 segfaults and 200 time-outs (signpost backend mostly). Ported tests
(`tests/reference_ports.rs`): `IntegerFlip`, `DelaunayTriangulationsAgree`, `IntegerTrace`,
`IntegerRefine`, `IntegerCommonSubdivision`, `CommonSubdivisionCompareIntegerSignpost`,
`FunctionTransfer`; not ported: the two trace-versus-bulk comparisons and
`CommonSubdivisionGeometry`.

## 11. Decisions

* **I1** shared connectivity (§2.1). **I2** `i64` with checked arithmetic; an overflow is
  `IntegerError::Overflow { edge }`, never a wrap (the paper notes counts can grow exponentially with
  Dehn twists, Sec. 3.7). **I3** local index start (§2.2). **I4** cost of the certificate, measured
  before choosing its default. **I5** clamped crossings. **I6** distance of a face split to the
  nearest arc. **I7** boundary removal. **I8** deterministic removal, refusal. **I9** tracing over
  `T¹` without stored signposts. **I10** lazy point location `T⁰ → T¹`. **I11** refinement order
  (as built: one shared driver and a shared removal choice, WP6). **I12** narrow-cone rule (decided:
  the reference's rule, I23; the reported minimum angle is over the faces the rule does not
  exclude, whereas the reference's statistic excludes faces with any narrow vertex). **I13** where
  mollification lives.

From the map of the reference (`runs_tmp/ic-reference-map.md`, "map §n"; `gc` paths abbreviated
`nc` = `src/surface/normal_coordinates.cpp`, `ic` = `src/surface/integer_coordinates_intrinsic_triangulation.cpp`,
`it` = `src/surface/intrinsic_triangulation.cpp`, `cs` = `src/surface/common_subdivision.cpp`):

* **I2, the reference.** `int` (32-bit), no overflow check, `int`/`size_t`/`double` mixed (map §1, §4).
  We keep `i64` checked.
* **I3, the reference.** Local indices CCW from geometry-central's `v.halfedge()`; at a boundary vertex
  the exterior halfedge gets `deg₀ − 1` (`nc:31-57`). Ours (I3) is the same order up to the start;
  the demo cannot dump `r`, so no parity is lost.
* **I14** *Eq. 1 as printed.* The reference's `cornerCoord` drops `max(0, ·)` and uses raw `n`, and its
  fan-versus-triforce test relies on a `size_t` promotion turning `/2` into a floor (`nc:790-792`,
  `nc:1094`; map §2.1). We compute Eq. 1 verbatim (`n⁺`, the `max`, an exact even numerator checked
  by `certify`) and decide "fan" by testing `e > 0` explicitly. Never `/` on a negative `i64`.
* **I15** *Face split counts are derived from a region, not from rounded slacks* (§4.4). The paper's
  Alg. 5/6 zero all slacks but the largest; the reference instead snaps the other corners to
  `c + σ_max` in the triforce case and uses a fixed precedence in the fan case, which can leave
  inconsistent counts when rounding puts the point in two corner bands (`nc:1107-1141`, map §5
  items 3–4). We enumerate the regions of the face's arc arrangement combinatorially (from the
  corner coordinates alone), use the line-side tests only to pick one region (the one whose test
  margins are largest, ties by region index), and compute `n_mi, n_mj, n_mk` from that region
  (`ic-face-split-normal`). Validity then holds by construction and is what Lean tier B proves.
* **I16** *App. A uses every corner of the region* (`3 ≤ ρ ≤ 6`), as the paper says; the reference
  picks 3–4 per case and rescues coordinates outside `[−0.5, 1.5]` by heuristics, then clamps
  (`ic:918-1158`). We solve the minimum-norm system (Eq. 15) by QR of the `3 × ρ` matrix, report the
  smallest barycentric coordinate before any clamp, and return `RecoverFailed` beyond `1e-12`.
  As built (WP4): an insertion place within `1e-12` outside `(0, 1)` is moved towards the region's
  centroid by the smallest fraction that makes it a valid `Face` point, and the fraction is
  reported (12 cases in the tests, smallest coordinate `−2.2e-13`); point queries instead snap
  coordinates `≤ 1e-12` to an edge or vertex point (`LOCATE_SNAP`, as §4.9). On flat meshes the
  recovered place is within `8.6e-14` mean edge lengths of the true one; the number of corners used
  always equals the region's (3/4/5/6).
* **I17** *Alg. 2 and 3 in full.* The reference omits Alg. 3's walk-back past inserted vertices
  (asserts the curve starts in `V⁰`, `ic:1980`) and Alg. 2's `uStart/uEnd` remapping of the input
  parameter for a piece of a split input edge (stores the piece's own parameter, `ic:292-299`; map
  §5 B4), and its continuation through inserted on-edge vertices only inspects the first corner
  (`ic:1924-1940`, B3). All latent in its pipeline (only boundary or marked edges are split). We
  implement Alg. 2–4 as printed, and the tests split interior input edges and then flip the
  pieces, which exercises exactly these paths.
* **I18** *Loops in tracing.* The reference reports a loop when the first and last crossing share a
  halfedge, without comparing indices (`nc:355-357`, B6), which misfires on spiralling curves. We
  detect a repeated `(halfedge, index)` pair, which cannot happen for a valid state (an error
  naming the crossing if it does).
* **I19** *No cached common subdivision.* The reference caches it and forgets to invalidate it after a
  face split (`ic:1184-1244`, B1). We compute it on request; the returned `CommonSubdivision`
  wraps the `Overlay` (which has no fingerprint field) with the state's fingerprint (FNV over
  faces, lengths, `n`, `r`, `q⁰`), so a stale one is detectable and `interpolation_matrices`
  refuses it.
* **I20** *Splits of non-shared fixed edges.* The reference throws after modifying the mesh when a
  marked edge with `n ≥ 0` is split (face split, then a refused flip; `ic:1479-1509`, B5). Our edge
  split with `n ≥ 0` performs the explicit edge split of the intrinsic edge (as for `n < 0`, with
  the crossings of `ij` divided at the split point by their strip parameters, which the reference
  has as dead code `nc:191-221`), so fixed edges never need the flip; a refused split changes
  nothing. Every edge split returns both halves (in order from the old edge's tail) and computes
  the lengths exactly as `ddg-intrinsic`'s `split_edge` does (IT-SPEC N10), so constrained paths
  split identically in both representations (FlipOut contract R4).
* **I21** *Self-edges.* The reference's `inEdge` and roundabout ordering assume distinct diamond
  vertices (B7). Every operation here is tested on Δ-complex configurations (self-edges, repeated
  diamond vertices) produced by flips on the generators' small meshes.
* **I22** *Mollification refreshes everything.* The reference mollifies after caching geometry (B8);
  here mollification happens before the triangulation is built (§7).
* **I23** *Refinement details of the reference, adopted.* The narrow-cone rule (I12) is the
  reference's, which is the paper's: skip a face with exactly one vertex of angle sum below `π/3`,
  or whose parent input face has such a vertex; skip degree-1 corners and corners between two
  fixed edges (`it:313-372`). The Dijkstra ball has radius `ℓ` of the split edge (`it:463-504`,
  Alg. 10 line 20). A circumcenter trace that fails skips the face, counted in the report
  (`trace_failures`); the reference skips silently (`it:538-542`). The order (I11) is the shared
  refinement driver's, i.e. IT-SPEC's, so the two representations of this workspace insert the same
  points; against the reference (max-heap by area, boundary first, up to 5 rescans) only counts
  and statistics are compared.
* **I24** *Removal choice* (I8 refined). The reference scores candidates by `min(A₁, A₂)/(A₁ + A₂)` of
  the two new triangles, prefers self-loops, gives up after `10 deg` flips, refuses boundary vertices
  (`ic:1603-1714`). We adopt the score and the self-loop preference (ties by edge index), the bound,
  and `RemovalStuck`; boundary inserted vertices: removal implemented (flip to one incident face,
  then merge the two boundary edges, both `n = −1` along the same input edge) only if a caller needs
  it, otherwise refused with `Unsupported` as in the reference (I7).
* **Departures of the reference kept for reference only**: mollification on and silent
  (`ic:38`; ours explicit, §7); `moveVertex` ignores a failed removal (ours refuses the move and
  restores); leftover debug output (`ic:773`) and `exit(1)` (`nc:930`) are not ported.
* **I25** *Face-split lengths.* Taken from `ddg_intrinsic::plan_face_split` (`it-point-distances`),
  not evaluated by Eq. 5, so that both representations insert the same points with bitwise the same
  lengths (I11). Eq. 5 (`ic-tangent-length`) stays a formula id: Lean proves it equal to
  `it-point-distances`.
* **I26** *Tangent frames for consumers (open, WP11).* `ddg-intrinsic-apps`' tangent fields
  (`vector/connection.rs`) use the signpost crate's stored input-aligned frames; `AngleFrames` are
  neither input-aligned nor stable under flips. For vertices of `V⁰` an input-aligned frame is
  available from the roundabouts plus one angle within the input wedge; the design is decided in
  WP11 before that consumer switches representation.
* **I27** *Removal of a vertex with a straight star (decided in WP3: `RemovalStuck`, state
  restored).* An interior vertex made by an edge split has degree 4 with two angles of `π`; those two
  spokes are refused by `plan_flip` as too thin, and the other two flip only when their
  quadrilateral is convex at the far end. Measured: right after the split 1728 of 2588 split
  vertices removed, 860 stuck; after Delaunay flipping 1692 removed, 1898 stuck (1873 straight
  stars). Refinement removes interior vertices, which are face insertions there, so it rarely meets
  this; FlipOut's Bézier splits and constrained paths do. **Future:** `Complex::merge_edge` (the
  inverse of `split_edge`) in `ddg-intrinsic`, so that such a vertex is removed by merging the two
  halves back.
* **I28** *Split point within rounding of a crossing (open).* The integer edge split refuses a
  point within rounding of a crossing (§4.5), so refinement skips such a face, where the signposts
  split and snap the place onto the input edge (1 of 256 differential runs, WP6). The alternative
  is to snap the integer split to the crossing (the new vertex on the input edge); undecided.
  The refine driver has no counter for these edge-split refusals (WP10); face-split refusals are
  `n_recover_failed` (7 on Thingi10k, in 2 meshes).
* **I29** *An integer simplicial export (open).* `ddg-uv`'s intrinsic path ends in
  `export_simplicial` and its cut transfer uses `SimplicialExport::edge_point` and
  `segments_on_input`; the integer representation has no export yet (it would be the common
  subdivision plus `extract_intrinsic_edge`), so `ddg-uv` stays on signposts until it exists.
* **I30** *Non-terminating refinement on Thingi10k 176193 and 190253 (open).* Both closed, about
  16k faces, with input angles near `1e-5°`. The shared driver `ddg_intrinsic::refine` keeps
  inserting: at 10 insertions per vertex 176193 stops after 162 790 insertions at minimum angle
  `6.7e-6°` with 18 671 corners below 25°; signposts make the identical insertions; the demo
  (integer backend) finishes them in 0.6 s and 1.4 s at 25.0°. The cause is in the driver, not the
  representation. Closing criterion: `tests/thingi_regressions.rs` (3 `#[ignore]`d tests, need
  `data/thingi10k/converted/`) pass.

## 12. Work packages

| WP | content | files | depends on | model |
|---|---|---|---|---|
| WP0 | **done** (2814e8a–aca460e): `Complex`, `Triangulation`, `plan_*`, `AngleFrames` (I9), `ddg_core::mollify` (I13); 240 bitwise pins unchanged | `crates/ddg-intrinsic/src/…`, `crates/ddg-core` | — | Opus |
| WP1 | **done** (e2d71da–2fb8762): `IntegerTriangulation` (no `Triangulation` impl yet: needs WP2–WP4), init, §2.3, the oracle's integer formulas, combinatorial Alg. 1–4, `certify`, `from_intrinsic` (conversion from signposts, the differential oracle for WP3), errors | `src/{lib,error,coords,formulas,trace,certify,convert}.rs` | WP0 | Opus |
| WP2 | **done** (7aa2891–bb56df2): geometric crossings, strip intersection, `extract_intrinsic_edge` (linear), `trace_edge_on_input` wired; same crossed input edges as signposts on all 678 806 edges, 0 clamps | `src/geometry.rs`, `tests/geometry.rs`, `examples/extract_cost.rs` | WP1 | Opus |
| WP3 | **done** (7b7dcd3–a315f9e): flip (Eq. 3, 4), Delaunay flipping, edge split incl. `n ≥ 0`, removal (I27), `impl Triangulation`; 437 595 flips on 266 states equal `from_intrinsic` of the signposts exactly | `src/{flip,split,remove,interface}.rs`, `tests/ops.rs` | WP1 | Opus |
| WP4 | **done** (548d804, 6080ce8): face split (I6, I15), App. A and the `T¹ → T⁰` point query (I16), moves; `insert_vertex` of the trait | `src/{insert,recover,move_vertex}.rs`, `tests/insert.rs` | WP2 | Opus |
| WP5 | **done** (432e0ea, 9cbd214): common subdivision, `P₀`, `P₁`, `T⁰ → T¹` point queries, batch `extract_all_edges` | `src/subdivision.rs`, `tests/subdivision.rs` | WP2 | Sonnet |
| WP6 | **done** (49bac6c, 91ca49b): Delaunay refinement with boundary, one driver over `Triangulation` (`RefineRules::{Ruppert, Chew}`) | `src/refine.rs`, `tests/refine.rs`, `crates/ddg-intrinsic/src/refine.rs` | WP3, WP4 | Opus |
| WP7 | **done** (b463c96, 406c0e0): transfer (Eq. 13–14); Fig. 14 as the paper's Poisson experiment (trends reproduced, factors depend on the unpublished generator) | `src/transfer.rs`, `tests/transfer.rs`, `examples/fig14.rs` | WP5 | Sonnet |
| WP8 | **done** (1149fdf–854e042; `ic-vertex-removal-delta` open, cited): Lean tiers A, B, B′; oracle | `lean/proofs/DdgProofs/IntegerCoords/…`, `lean/oracle/DdgOracle/IntegerCoords…`, `CITED.md` section | this file | Opus |
| WP9 | Lean tier C, the complex model | `lean/proofs/DdgProofs/IntegerCoords/Complex/…` | WP8 | Opus |
| WP10 | **done** (307308a–a9fa390): Thingi10k (Table 1 integers 99.97 %, two time-outs, I30), benches, the demo comparison, ported reference tests (§10.2, §10.3) | `examples/{thingi_convert,thingi_eval}.rs`, `benches/representation.rs`, `tests/{reference_ports,thingi_regressions}.rs`, `scripts/{ic-eval,ic-demo-eval,ic-demo-build,ic-timings}.sh` | WP6, Q1, Q2 | Sonnet |
| WP11 | **done, default still signposts** (b1fe875–bebb60b): `ddg_intrinsic::generic_ops` (`cotan_laplacian`, `lumped_mass`, bitwise as the signpost methods); `ddg-intrinsic-apps` fem, steiner, amr, heat generic (`Representation`, `InputPoints`); Python `IntrinsicTriangulation(…, representation="integer")` with refine, common subdivision, edge polylines, transfer, certify, and `representation=` on `geodesic_distance`, `steiner_tree`, `adaptive_solve`; Blender "Show Intrinsic Refinement" (`intrinsic_draw.py` reusable). Still signposts-only: `optimal_delaunay`, `move_vertex`, `odt_target`, `export_simplicial`, `fingerprint`, `SimplicialExport::edge_point`, `segments_on_input` (so `ddg-uv`, I29), `signpost` (tangent fields, I26) | per crate | WP5, WP6 | Sonnet |

At most three agents at a time; timings alone.

## 13. The user's decisions (2026-10-31)

* **Q1** Thingi10k: download all of it. Raw archive `Thingi10K.tar.gz` (9.6 GB) and the metadata CSVs
  at Hugging Face revision `v1.5.0` into `data/thingi10k/` (gitignored). The paper converted the raw
  files with MeshLab to PLY; our conversion (WP10, `examples/thingi_convert.rs`) merges exactly
  coincident vertices, keeps the meshes that are manifold triangle meshes after that, and lists
  every excluded file with the reason (`data/thingi10k/converted.tsv`). As built: 7745 kept
  (paper 7696; the 49 cannot be attributed without MeshLab's rules); excluded: non-manifold
  vertex 649, duplicate face 616, repeated vertex in a face after the merge 482, non-manifold edge
  343, inconsistent orientation 160 (95 repairable, not repaired), no faces 3, truncated ASCII STL
  2. Of the metadata's 7754 manifold meshes 7744 are kept; one kept mesh (73178) is non-manifold
  in the metadata. Downloaded and extracted 2026-10-31
  (`data/thingi10k/raw/Thingi10K/raw_meshes/`, 10 000 files, 37 GB): the dataset's own metadata
  (`geometry_data.csv`, 9 997 rows) marks 7 754 meshes vertex- and edge-manifold (all oriented;
  7 346 closed).
* **Q2** Build the demo: yes. Polyscope (`56b82cd`, with its submodules) and happly (`88f7972`)
  fetched into the demo's `deps/`; geometry-central at `ab116f1`. Built outside the repository,
  `--noGUI` (WP10). Recipe `scripts/ic-demo-build.sh`, build in `runs_tmp/ic-demo-build/`:
  `GC_EIGEN_LOCATION`, `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` (polyscope's glad), the demo's
  `deps/happly` on `CPLUS_INCLUDE_PATH` (geometry-central's copy is an empty submodule).
* **Q3** Consumers use integer coordinates by default once the evaluation passes; signposts stay an
  option (WP11). As built: the switch is `Representation::DEFAULT` in
  `crates/ddg-intrinsic-apps/src/representation.rs` plus the `representation=` defaults of the
  `ddg-py` bindings (still `Signposts`); flipping it after WP10 is one commit. Checklist:
  `Representation::DEFAULT`; the literal `representation="signposts"` defaults of
  `geodesic_distance`, `steiner_tree`, `adaptive_solve` (FlipOut's bindings follow `DEFAULT`
  through `None`); Blender `geodesic_operator.REPRESENTATION_DEFAULT`; stubs and READMEs; the tests
  that assert the default (`the_default_representation_is_signposts`, `test_integer_coords.py`,
  `test_default_representation_is_the_documented_one`). `ddg_tools.IntrinsicTriangulation` and
  "Show Intrinsic Refinement" already default to integers, by design (they exist to show them).
* **Q4** FlipOut (Sharp and Crane 2020a) is in scope, as its own crate `ddg-flip-geodesics` with its
  contract `crates/ddg-flip-geodesics/SPEC.md`, on top of the shared triangulation interface
  (PDF downloaded to `../FlipOut/FlipOut.pdf`; reference `gc: src/surface/flip_geodesics.cpp`).
* **Q5** Blender gets a visible operator: the intrinsic Delaunay refinement / common subdivision as a
  mesh object (Fig. 12 style), besides the internal switch of Q3 (WP11). The FlipOut operators
  reuse its drawing of intrinsic edges on the input mesh.
* **Q6** Lean tier C: yes ("include the Lean combinatorial model … prefer to have the relevant
  theorems be proved"); WP9, after WP8. A survey of existing formalizations of every theorem of §9
  (`runs_tmp/ic-formalization-survey.md`) decides how tiers C–E are built.
