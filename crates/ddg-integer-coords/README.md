# ddg-integer-coords

Intrinsic triangulations stored by integer coordinates, after

> Mark Gillespie, Nicholas Sharp, Keenan Crane. *Integer Coordinates for Intrinsic Geometry
> Processing.* ACM Transactions on Graphics 40(6), Article 252, 2021.
> <https://doi.org/10.1145/3478513.3480522>

The contract is [`SPEC.md`](SPEC.md) (decisions `I1`–`I27`, work packages in §12). This crate is
the second representation of intrinsic triangulations in the workspace; the first, signposts
(Sharp, Soliman, Crane 2019), is `crates/ddg-intrinsic`, whose `Complex` (the connectivity) it
shares.

## What is here (WP1)

| | |
|---|---|
| `IntegerTriangulation` | connectivity (`ddg_intrinsic::Complex`), `Arc<InputSurface>`, lengths `f64`, normal coordinates `n: i64` per edge (SPEC I2), roundabouts `u32` per halfedge (`NO_ROUNDABOUT` for tails in `V*`), places `q⁰` of inserted vertices. `new(mesh, geometry)` is `T¹ = T⁰`; `with_options(mesh, lengths, &InitOptions)` mollifies first if asked (`ddg_core::mollify`, off by default, `InitReport::mollified`). `from_parts` / `into_parts` for stored data. Read-only accessors. |
| `LocalIndices` | the local indices of Supp. B with the start of SPEC I3 (`TriMesh::outgoing_halfedge`; at a boundary vertex the last index is the arriving boundary edge, no wrap). |
| `formulas` | Eq. 1 (`corner_coordinate`, never a floor division, I14), Eq. 2 (`emanating`), face validity, `reverse_crossing`, Alg. 1's step (`trace_step`), Alg. 3's index (`identify_index`), Eq. 6 (`split_edge_normal`), `FaceCoordinates`; all equal to the Lean oracle's. |
| `trace` | combinatorial Supp. Alg. 1–4: `trace_from` / `trace_with` (Alg. 1), `extract_curve` (Alg. 2 without geometry), `identify_traced_edge` (Alg. 3, walking back past inserted vertices on the edge, halfedge `a→i`, SPEC §3.1), `extract_edge` (Alg. 4 with the continuation through inserted vertices, I17), `input_directions` (the window of Alg. 4 at a vertex, all directions at once), `curve_ends`. Generic over `NormalGluing` so that hand-built Δ-complexes can be traced. |
| `certify()` | SPEC §2.4, all four items (below). |
| `from_intrinsic` | integer coordinates read off a signpost triangulation by tracing its edges (below). |

**Not yet** (later work packages): face split and App. A (WP4; `insert_vertex` of the shared
interface is `Unsupported`), the common subdivision (WP5), refinement (WP6), transfer (WP7).

## Operations (WP3)

| | |
|---|---|
| `flip_edge`, `can_flip` | refusals and new length from `ddg_intrinsic::plan_flip` (bitwise as signposts); Eq. 3 with the corrected first terms `c^{ji}_l + c^{ij}_k` in its doubled form (`i128`, an odd value is an error, I14; outside `i64` → `Overflow { edge }`; `n_kl < −1` → `InvalidCoordinates` naming the edge); Eq. 4 for both new halfedges with tails in `V⁰` (`LocalIndices::modulus`, no reduction at the boundary). All from the pre-flip state, then `Complex::flip_edge` + `move_halfedge_data`. |
| Delaunay flipping | the shared `Triangulation::flip_to_delaunay`; `flip_to_delaunay_certified` certifies after it (I4). |
| `split_edge` | refusals and lengths from `plan_edge_split`. `n = −1` (all boundary edges): Eq. 6 on both sides, `q⁰ = (1 − t) q⁰_i + t q⁰_j` on the input edge (the signposts' arithmetic). `n ≥ 0` (I20, explicit, also for constrained edges): crossings before the point from `extract_intrinsic_edge` (`distance < t ℓ`), `q⁰` interpolated in the input face between the two neighbouring crossings (refused within rounding of a crossing), the integers by `formulas::split_side` (the arcs separated by the new chord; Eq. 6 is its `n = −1` case). `split_edge_given` takes both decisions from the caller. |
| `remove_vertex` | inserted interior vertices: flips by I24 (degree must drop; loops first, then `min(A₁, A₂)/(A₁ + A₂)`, ties by edge id; bound `10 deg`), then `Complex::remove_vertex` with `VertexRemoval` moving `n`, `r`, lengths, `q⁰`. Input vertex / constrained spoke: refused before anything changes; boundary: `Unsupported` (I7); no flip left: `RemovalStuck`, with the state restored exactly (connectivity from a copy, the rest from a journal). Returns `Removal` (the flips, the renumbering). |
| `impl Triangulation` | lengths, angle sums, places, flips, splits, the constrained mask; `trace_edge_on_input` = `extract_intrinsic_edge` (WP2); errors through `to_intrinsic_error` (`IntrinsicError::Representation`, `"integer coordinates"`). |

**I27.** A vertex made by an edge split has two spokes at an angle of `π`: their flips would make
a triangle without area (`plan_flip` refuses them), and of the two halves the one with the angle
below `π` at the vertex can be flipped only if its quadrilateral is convex at the far end. When
it cannot, the removal is `RemovalStuck` (state restored); merging the halves would need an
inverse of `Complex::split_edge`, which `ddg_intrinsic` does not have. Measured below: about a
third of fresh split vertices, half after Delaunay flipping; face insertions almost never (12 of
1 160).

**Timings** (release, one pinned core, load average 3.3, 2026-10-31, best of 5 / 3;
`cargo test --release --test ops -- --ignored --nocapture timing`):

| mesh | \|E\| | Delaunay flipping (flips) | signposts | integer | certify after | random flips/s, signposts | integer | removal |
|---|---|---|---|---|---|---|---|---|
| spot | 8 784 | 281 | 1.57 ms | 1.63 ms | 1.48 ms | 3.74 M | 3.83 M | 7.8 µs |
| camelbad | 34 116 | 8 110 | 8.87 ms | 9.29 ms | 6.94 ms | 3.27 M | 3.23 M | 25 µs |

Integer flips cost the same as signpost flips (Eq. 3, 4 are a few dozen integer operations; the
time is in `plan_flip`, shared). A removal copies the connectivity once (O(|E|), most of its
cost on camelbad) so that a stuck removal can be undone exactly.

## Geometry of the crossings (WP2)

Connectivity comes from the integers; floating point decides only where the crossings lie. Both
directions unfold a strip of triangles (`it-unfold-apex`) and cut one straight segment with its
interior edges (`it-ray-edge-intersection`); a parameter outside `[0, 1]` by rounding is clamped
and counted (`StripReport::clamped_crossings`, SPEC I5).

| | |
|---|---|
| `extract_geometric_curve(ζ)` | Supp. Alg. 2: the segment through a crossing, its input edge (Alg. 3), `u_start`/`u_end` from the places of its ends (I17), `GeometricCrossing { crossing, input_edge, u, v, sin_angle }` (`u` in the input edge's direction, `v` along the crossed halfedge, `sin_angle` the conditioning). `extract_geometric_crossing(ζ)`: one crossing of it. `segment_geometry`: the strip of a given `CurveSegment`. |
| `extract_edge_geometry(h⁰)` | Alg. 4 with positions: WP1's `InputEdgePath`, the parameter of every vertex on it, the crossings of every piece. |
| `extract_intrinsic_edge(e, &EdgeTraceOptions)` | `ExtractIntrinsicEdge` (FlipOut R2), returned as the signposts' `EdgeTrace`: each crossing's input edge identified from the integers (Alg. 1 + 3, memoized along the walks), the input faces followed exactly from the tail's roundabout or `q⁰` to the head's (consecutive input edges must bound a common face), positions from the straight segment through those input faces. `mismatch` and `clamp` are reported, never a failure. `_with_report` adds the `StripReport`. It is `trace_edge_on_input` of the shared interface. |

Tests (`tests/geometry.rs`, the 266 states of `tests/states.rs`, 37 of them with self-edges):
`extract_intrinsic_edge` of all 678 806 edges crosses exactly the input halfedges the signposts'
`trace_edge` crosses, in the same order (272 651 crossings, first faces and `along` equal); the
two directions of every crossing (Alg. 2 over `T¹`, `extract_intrinsic_edge` over `T⁰`) and the
signposts' overlay agree to `|Δ| sin θ ≤ 3.5e-10`; raw `|Δt|` exceeds `1e-9` at 2 grazing
crossings only (in 3 states; sin θ = 1.4e-4 and 9.3e-4, `|Δt|` up to 5.1e-8); 0 clamps; mutated integers give
results or errors, never panics.

Cost (`examples/extract_cost.rs`, after Delaunay flipping, release, one pinned core, load
average 10.7 from other jobs, best of five): all edges by `extract_intrinsic_edge` spot 0.33 ms
(37 ns per edge + crossing), camelbad 5.3 ms (112 ns); the signposts' `trace_edges` 0.20 / 3.3 ms.
All input edges by `extract_edge_geometry` 7.0 / 31 ms, of which WP1's combinatorial
`extract_edge` is 6.4 / 27 ms (`input_directions` at both ends of every edge).

## The certificate

`certify()` returns a `CertifyReport` or `IntegerError::InvalidCoordinates` with the item of
SPEC §2.4 (`ValidityItem`), the place (face, edge, halfedge, vertex, crossing or input edge) and
the values. In order:

0. ranges: array sizes, `n ≥ −1`, `n = −1` on boundary edges, roundabouts in `[0, deg₀)` exactly
   on halfedges leaving input vertices, `q⁰` of the right kind, boundary status of input vertices;
1. item 1: even numerators of Eq. 1 in every face (item 2 follows from Eq. 2);
2. item 3: every edge with `n = −1` joins two distinct vertices along one input edge, named by the
   roundabouts at its input-vertex ends and by the places of its inserted ends;
3. item 4: at every input vertex the roundabouts agree with the emanating counts (going ccw, the
   next roundabout is `r + n⁻ + e`, the counts sum to `deg₀`); every input edge followed from its
   tail arrives at its head by the end the roundabouts there give, through inserted vertices that
   lie on it in order (two curve ends each); every crossing and every `n = −1` edge is used by
   exactly one input edge (a bitmap of `Σ n⁺` bits).

**Cost (SPEC I4).** Linear: every crossing is visited once. `examples/certify_cost.rs`, release,
one pinned core, load average 0.6 (2026-10-31), best of five, after `k·|E|` random flips:

| mesh | flips | \|E\| | Σn⁺ | certify | ns / (\|E\| + Σn⁺) | ≈ flips | signposts' `trace_edges` |
|---|---|---|---|---|---|---|---|
| icosphere 4 | 0 | 7 680 | 0 | 1.2 ms | 159 | – | 0.1 ms |
| icosphere 4 | 50 E | 7 680 | 11 704 | 2.0 ms | 101 | 13 400 | 2.6 ms |
| torus 48×16 | 50 E | 2 304 | 4 708 | 0.59 ms | 84 | 4 300 | 0.8 ms |
| bumpy grid 40×40 | 50 E | 4 880 | 6 477 | 1.1 ms | 101 | 8 700 | 1.5 ms |
| spot | 10 E | 8 784 | 9 876 | 2.1 ms | 114 | 14 600 | 2.5 ms |
| spot | 50 E | 8 784 | 13 353 | 2.3 ms | 103 | 15 800 | 3.0 ms |
| camelbad | 50 E | 34 116 | 55 244 | 9.8 ms | 110 | 60 000 | 13.1 ms |

"≈ flips" is the certificate's time over the time of one signpost flip on the same mesh (0.13 to
0.16 µs): one certificate costs 1.7 to 1.9 `|E|` flips, a little less than tracing every edge with
signposts. Recommended default: certify at the end
of every public operation that does a variable amount of work (Delaunay flipping, refinement,
before the common subdivision, a FlipOut path), never per flip or split; in debug builds and tests,
after every operation on small meshes and per batch on large ones.

## Conversion from signposts

`IntegerTriangulation::from_intrinsic(&IntrinsicTriangulation, &EdgeTraceOptions)` traces every
edge over the input mesh with the signposts: `n` is the number of input edges crossed (`−1` along
an input edge), the roundabout of an end at an input vertex is the local index after the input
corner where the path leaves it, `q⁰` is the signposts' `vertex_location`. Floating point decides
this, so the result is certified before it is returned; a trace that fails or integers that do not
certify are an error. Two exact readings are applied before: a path without crossings between two
points of one side of its input face runs along that side, and a last crossing with the input edge
on which the far vertex lies is the arrival. On the 17 standard shapes of `ddg-intrinsic`'s tests,
spot and camelbad, in 7 scenarios (initial, random flips, Delaunay, insertions, edge splits,
splits of input edges then flips, mixed) and 2 seeds, **266 of 266 conversions certified**; neither
exact reading was needed.

## Tests

* `tests/states.rs`: the 266 conversions; every crossing of every edge (272 651), extracted in both
  directions and identified by Alg. 3, names the input edge the signposts' trace crosses there
  (1 328 of them by walking back past inserted vertices); `extract_edge` of every input halfedge
  arrives where it should and the paths use `Σ n⁺` crossings; **15 941 single mutations** (`n ± 1`,
  a roundabout `± 1`, two roundabouts at a vertex swapped) all rejected — by item 4 (extraction):
  6 167, item 3: 7 414, ranges: 2 240, item 1: 120. Run with `--nocapture` for the tables.
* `tests/oracle.rs`: the Lean oracle's committed fixtures (`lean/oracle/fixtures/ic`), exact
  equality: `corner` (519 faces), `trace` (192 steps), `tracecurve` (14 curves on the Fig. 5 complex
  and two hand-built Δ-complexes, 2 closed curves), `identify` (168 indices; 34 arguments of the
  random fixture have an invalid face `ija`, where the oracle rounds Eq. 1 and this crate refuses,
  I14), `splitedge` (52), `flip` (305 quads, 302 roundabout updates, 304 flips back).
* `tests/ops.rs` (WP3), against the signposts on the 17 shapes, spot and camelbad, 7 scenarios, 2
  seeds; after every batch the integer state **equals** `from_intrinsic` of the flipped signposts
  (connectivity, lengths bitwise, every `n`, roundabout and place) and certifies:
  * flips: 266 states, 1 008 batches, 912 744 attempts, 437 595 flips, 475 149 refusals (all the
    same in both), then Delaunay flipping with the identical flip sequence and report (183 002
    flips); Δ-complexes up to 23 self-edges and 46 faces with a repeated vertex in one state
    (I21); 0 mismatches;
  * flip then flip back: 10 490 edges, integers restored exactly;
  * edge splits: 4 256 along input edges known to the signposts (exact, places bitwise), 1 060
    along input edges the signposts trace (all but the place exact, places within 2e-11), 5 030
    across input edges with the signposts' crossings and place (exact) and the same 5 030 with this
    crate's own crossings (all but the place exact, places within 1.5e-10 barycentric); refusals
    identical;
  * removals of inserted vertices, the same flips and removal replayed in the signposts, compared
    exactly: 3 884 removed (9 416 flips, 2.4 per removal, at most 14), 1 301 stuck (1 287 of them at
    a straight star from an edge split), states restored exactly; refused removals (input vertex,
    boundary, constrained spoke) change nothing;
  * I27: 2 588 vertices from splits of interior edges (913 across input edges) removed at once:
    1 728 removed, 860 stuck; after Delaunay flipping 1 692 removed, 1 898 stuck (1 873 straight);
  * overflow: hand-built states near `i64::MAX` give `Overflow`, never wrap; Eq. 3 below `−1`
    (quad `2,0,0,0,0`) is refused as invalid.

## The reference implementation

geometry-central (MIT, © 2017–2019 Nicholas Sharp and the geometry-central contributors), commit
`ab116f1`, pinned by the authors' demo, was read for the details the paper leaves out (map:
`runs_tmp/ic-reference-map.md`); nothing was translated from it. Its notice is in
`LICENSES/geometry-central-MIT.txt`.
