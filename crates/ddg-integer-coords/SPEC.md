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
subdivision, a FlipOut run) when `IntegerOptions::certify` is set (default on); never per flip or
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
state, then `complex.flip_edge(e)`, whose `ComplexFlip` says where `n_kl`, `r_kl`, `r_lk` go. Exact: the
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
statistics leave out.

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

### 4.5 Edge split (Sec. 3.4, Alg. 7, Eq. 6)

`n_ij ≥ 0`: face split then flip (Fig. 10 top). `n_ij = −1`: explicit split along the shared
curve, `n_mj = n_mi = −1`, `n_mk = max(n_ki, n_jk, 0)` (`ic-split-edge-normal`; the `0` matters when both are `−1`), and the same on the
other side when interior; `q⁰_m = u_i q⁰_i + u_j q⁰_j` along the input edge. Boundary edges are
always of the second kind. Lengths and refusal from `plan_edge_split`, connectivity from
`Complex::split_edge` (I20: also for `n_ij ≥ 0`, without the face-split-then-flip route).

### 4.6 Vertex removal (Sec. 3.5, Alg. 8, Thm. D.1)

Only `i ∈ V*`. Flip flippable edges at `i` until its degree is 3 (interior) or the boundary
analogue (I7: the reference's rule for boundary vertices, `runs_tmp/ic-reference-map.md`), then
delete. No integer data changes beyond the flips. Thm. D.1 proves that a simplicial star with
angle sum `2π` and degree `≥ 4` has a spoke whose quadrilateral has all angles `≤ π`; it need not
be strictly convex (Lean `cross_star_no_strict`: in the degree-4 cross every inner corner is `π`),
so N5's strict flippability can fail, e.g. at an interior vertex made by an edge split (**I27**).
For a star that is not simplicial the paper conjectures it (open, `ic-vertex-removal-delta`).
**I8:** the flip choice is deterministic (the flippable edge that lowers the degree with the
largest minimum angle afterwards; ties by edge index), the loop is bounded by `deg(i)²` flips, and
a star where no edge is flippable is refused with `RemovalStuck { vertex, degree }`, the state
restored. Counted over every removal of the evaluation (§10).

### 4.7 Moving an inserted vertex (Sec. 3.6)

`x = exp_i(v)` traced on `T¹` (IT-SPEC tracing over the intrinsic triangulation needs signposts;
here **I9:** the tangent vector is laid out in the faces around `i` from the edge lengths, i.e.
the trace is `ddg_intrinsic::AngleFrames::of(&t).trace(start, dir, len, &opts)`, the same kernel
with signposts computed from corner angles. Frame convention: direction 0 at a vertex is
`Complex::outgoing_halfedge`, the others the running sum of corner angles times `2π/Θ` (`π/Θ` at
the boundary), `Θ` from the lengths; a frame is valid only for the current connectivity, since a
flip at an interior vertex can move its outgoing halfedge. Measured against stored signposts:
within 6.2e-12 rad on 17 shapes, traces end within 6.2e-12 of the mean edge length), then `SplitFace(x)`, then `RemoveVertex(i)`. The optimal Delaunay
iteration of IT-SPEC N21 then runs on this representation unchanged.

### 4.8 Common subdivision (Sec. 3.7, Supp. Alg. 9)

Per face of `T¹`, oriented so that `n_ij ≥ n_jk, n_ki` (Alg. 9 line 5): case 1 (`e^{ij}_k = 0`,
corner-to-corner pairing) or case 2 (fan from the one corner with `e > 0`), Fig. 11. Vertex
positions: `V¹` from `q⁰` and the vertex itself, crossings from `ExtractEdge` over every `ab ∈ E⁰`
(Alg. 9 lines 15–20). Output: `ddg_intrinsic::Overlay` (the same type, so IT-SPEC N12's tests and
every consumer apply unchanged), plus the interpolation matrices `P₀` (`|V^S| × |V⁰|`) and `P₁`
(`|V^S| × |V¹|`), rows of 1 to 3 entries summing to one (Sec. 4.3). A face whose arcs fail the
pairing is impossible for a valid state; if it happens it is an error naming the face (it would
mean §2.4 was violated, which `certify` must have caught).

### 4.9 Point queries

`T¹ → T⁰`: a point of a face of `T¹` → the region of the face → `RecoverBarycentric` (as §4.4
without inserting). `T⁰ → T¹`: by the common subdivision of the containing input face (I10: lazily,
per input face, from `ExtractEdge` of its three edges). Same signatures as `ddg-intrinsic`'s
`point_on_input` / `point_on_triangulation`.

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
  evaluated without boundary.
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
* **Fig. 14** reproduced exactly as described (§13 has no question: the generator is ours, the
  random seeds are listed).
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

## 11. Decisions

* **I1** shared connectivity (§2.1). **I2** `i64` with checked arithmetic; an overflow is
  `IntegerError::Overflow { edge }`, never a wrap (the paper notes counts can grow exponentially with
  Dehn twists, Sec. 3.7). **I3** local index start (§2.2). **I4** cost of the certificate, measured
  before choosing its default. **I5** clamped crossings. **I6** distance of a face split to the
  nearest arc. **I7** boundary removal. **I8** deterministic removal, refusal. **I9** tracing over
  `T¹` without stored signposts. **I10** lazy point location `T⁰ → T¹`. **I11** refinement order.
  **I12** narrow-cone rule. **I13** where mollification lives.

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
  smallest barycentric coordinate before any clamp, clamp only coordinates within `1e-12` of 0, and
  return `RecoverFailed { face }` otherwise.
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
  face split (`ic:1184-1244`, B1). We compute it on request; a returned `Overlay` carries the
  state's fingerprint (IT-SPEC N15) so a stale one is detectable.
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
* **I27** *Removal of a vertex with a straight star (open, WP3).* An interior vertex made by an
  edge split has degree 4 with two angles of `π`, so no spoke is strictly flippable (§4.6). WP3
  decides between undoing the split (the inverse of `Complex::split_edge`, merging the halves) and
  `RemovalStuck`, and counts how often refinement meets it.

## 12. Work packages

| WP | content | files | depends on | model |
|---|---|---|---|---|
| WP0 | **done** (2814e8a–aca460e): `Complex`, `Triangulation`, `plan_*`, `AngleFrames` (I9), `ddg_core::mollify` (I13); 240 bitwise pins unchanged | `crates/ddg-intrinsic/src/…`, `crates/ddg-core` | — | Opus |
| WP1 | **done** (e2d71da–2fb8762): `IntegerTriangulation` (no `Triangulation` impl yet: needs WP2–WP4), init, §2.3, the oracle's integer formulas, combinatorial Alg. 1–4, `certify`, `from_intrinsic` (conversion from signposts, the differential oracle for WP3), errors | `src/{lib,error,coords,formulas,trace,certify,convert}.rs` | WP0 | Opus |
| WP2 | **done** (7aa2891–bb56df2): geometric crossings, strip intersection, `extract_intrinsic_edge` (linear), `trace_edge_on_input` wired; same crossed input edges as signposts on all 678 806 edges, 0 clamps | `src/geometry.rs`, `tests/geometry.rs`, `examples/extract_cost.rs` | WP1 | Opus |
| WP3 | flip (Eq. 3, 4 in `formulas.rs`, checked or `i128`, against the oracle's `flip` fixtures), Delaunay flipping, edge split, removal (I27); differential test: flip the same edges in both representations and compare with `from_intrinsic` exactly | `src/{flip,split,remove}.rs`, `tests/ops.rs` | WP1 | Opus |
| WP4 | face split, App. A, moves | `src/{insert,recover,move_vertex}.rs`, `tests/insert.rs` | WP2 | Opus |
| WP5 | common subdivision, `P₀`, `P₁`, point queries | `src/subdivision.rs`, `tests/subdivision.rs` | WP2 | Sonnet |
| WP6 | Delaunay refinement with boundary | `src/refine.rs`, `tests/refine.rs` | WP3, WP4 | Opus |
| WP7 | transfer (Eq. 13–14), Fig. 14 | `src/transfer.rs`, `examples/fig14.rs` | WP5 | Sonnet |
| WP8 | **done** (1149fdf–854e042; `ic-vertex-removal-delta` open, cited): Lean tiers A, B, B′; oracle | `lean/proofs/DdgProofs/IntegerCoords/…`, `lean/oracle/DdgOracle/IntegerCoords…`, `CITED.md` section | this file | Opus |
| WP9 | Lean tier C, the complex model | `lean/proofs/DdgProofs/IntegerCoords/Complex/…` | WP8 | Opus |
| WP10 | evaluation: Thingi10k, benches, the reference comparison | `examples/thingi.rs`, `benches/`, `scripts/ic-eval.sh` | WP6, Q1, Q2 | Sonnet |
| WP11 | consumers: `ddg-intrinsic-apps`, `ddg-uv`, Python, Blender (§13 Q3, Q5). Not in the trait yet (WP0 survey): `cotan_laplacian`, `lumped_mass` (could be provided methods); `refine`, `refine_with`; `optimal_delaunay`, `move_vertex`, `odt_target`; `export_simplicial`, `fingerprint`, `SimplicialExport::edge_point`, `segments_on_input`; `point_on_input`, `point_on_triangulation`; `overlay`, `check`; `signpost` (I26) | per crate | WP5, WP6 | Sonnet |

At most three agents at a time; timings alone.

## 13. The user's decisions (2026-10-31)

* **Q1** Thingi10k: download all of it. Raw archive `Thingi10K.tar.gz` (9.6 GB) and the metadata CSVs
  at Hugging Face revision `v1.5.0` into `data/thingi10k/` (gitignored). The paper converted the raw
  files with MeshLab to PLY; our conversion (WP10) merges exactly coincident vertices, keeps the
  meshes that are manifold triangle meshes after that, and lists every excluded file with the
  reason; the count is compared with the paper's 7696. Downloaded and extracted 2026-10-31
  (`data/thingi10k/raw/Thingi10K/raw_meshes/`, 10 000 files, 37 GB): the dataset's own metadata
  (`geometry_data.csv`, 9 997 rows) marks 7 754 meshes vertex- and edge-manifold (all oriented;
  7 346 closed).
* **Q2** Build the demo: yes. Polyscope (`56b82cd`, with its submodules) and happly (`88f7972`)
  fetched into the demo's `deps/`; geometry-central at `ab116f1`. Built outside the repository,
  `--noGUI` (WP10).
* **Q3** Consumers use integer coordinates by default once the evaluation passes; signposts stay an
  option (WP11).
* **Q4** FlipOut (Sharp and Crane 2020a) is in scope, as its own crate `ddg-flip-geodesics` with its
  contract `crates/ddg-flip-geodesics/SPEC.md`, on top of the shared triangulation interface
  (PDF downloaded to `../FlipOut/FlipOut.pdf`; reference `gc: src/surface/flip_geodesics.cpp`).
* **Q5** Blender gets a visible operator: the intrinsic Delaunay refinement / common subdivision as a
  mesh object (Fig. 12 style), besides the internal switch of Q3 (WP11). The FlipOut operators
  reuse its drawing of intrinsic edges on the input mesh.
* **Q6** Lean tier C: yes ("include the Lean combinatorial model … prefer to have the relevant
  theorems be proved"); WP9, after WP8. A survey of existing formalizations of every theorem of §9
  (`runs_tmp/ic-formalization-survey.md`) decides how tiers C–E are built.
