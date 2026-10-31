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

**Not yet** (later work packages): `impl Triangulation for IntegerTriangulation` (needs flips,
splits and extraction of positions), flips and Delaunay flipping (WP3), edge split and removal
(WP3), face split and App. A (WP4), the geometry of the crossings and strip unfolding (WP2),
the common subdivision (WP5), refinement (WP6), transfer (WP7).

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
  I14), `splitedge` (52).

## The reference implementation

geometry-central (MIT, © 2017–2019 Nicholas Sharp and the geometry-central contributors), commit
`ab116f1`, pinned by the authors' demo, was read for the details the paper leaves out (map:
`runs_tmp/ic-reference-map.md`); nothing was translated from it. Its notice is in
`LICENSES/geometry-central-MIT.txt`.
