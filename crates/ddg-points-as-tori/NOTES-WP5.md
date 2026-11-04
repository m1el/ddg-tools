# WP5 notes: baselines and mesh ground truth

Branch `points-as-tori`, build dir `target-pat5`, logs in `runs_tmp/pat-wp5/`. No GPU used.

## What exists

| file | content |
|---|---|
| `src/baselines/areas.rs` | Barill point areas: `point_areas`, `PointAreaOptions { k: 20, filter_opposite: true, threads }`, `PointAreaReport` (clipped, filtered, coincident, isolated, total), `voronoi_cell_area` `[pat-point-area]` |
| `src/baselines/winding.rs` | point-cloud GWN by direct sum (`winding_number(s)`, `p = x` terms left out and counted) `[pat-winding-dipole]`; regularized WN (Chen et al. 2024) with `S(t) = P(3/2, t²)` by series / Legendre continued fraction (no cancellation at small `t`, finite at `ρ = 0`) `[pat-regularized-wn]`; `winding_signed_distance` |
| `src/baselines/fwn.rs` | Barill FWN: octree for points (leaf 2, as the paper), median-split AABB tree for triangles (leaf 8), order 0/1/2 Taylor expansion, accuracy `β` (default 2) `[pat-fwn-expansion]`; triangle moments exact by the edge-midpoint rule; leaves exact (dipoles / van Oosterom–Strackee solid angles) |
| `src/baselines/kernel.rs` | `KernelField`: SHC (Eq. 28), SSPD (Eq. 29), LSE on the PAT neighbor sets and λ rules (reuses `blend::shift`, `Neighborhood`, `LambdaRule`, `mean_neighbor_distance`), frozen gradients, batch on the explicit pool `[pat-shc, pat-shc-gradient, pat-sspd, pat-sspd-gradient, pat-lse]` |
| `src/truth.rs` | `MeshTruth`: exact distance by an AABB tree + Ericson closest point (degenerate triangles by their edges; ties to the smaller face), sign by exact GWN (≤ 20k faces) or FWN (`SignMode::{Exact, Fast{β}, Auto}`), threshold ½ or `JAX_THRESHOLD = ¼`; `check_mesh` (boundary / non-manifold / inconsistent edges, degenerate faces, signed volume); `require_closed` (default) rejects open, non-manifold, inconsistent or inward meshes, naming the edge and faces |
| `scripts/pat_spsr.py` | SPSR via Open3D (`uv run` with PEP 723 deps, Python 3.12, open3d 0.20.0): depth = ⌈log₂ grid⌉ (6 for 64³), width 0, scale 1.1, single thread (paper), no trimming; writes `vertices.npy`, `faces.npy`, `densities.npy`, `meta.json` |
| `tests/baselines.rs` | 17 tests + the ignored `sanity_grid_errors` |

## Tests (`cargo test --release -p ddg-points-as-tori --test baselines`, 0.3 s; all pass, clippy clean for these files)

- Voronoi cell: hexagon exact to 1e-14; one-sided neighbors clipped.
- Point areas: plane grid interior cells `h²` to 1e-14 (only boundary cells clipped); Fibonacci sphere Σ A within 2 % (2000 pts) / 1 % (8000), no clipping; duplicates share the cell (½ each); the opposite sheet is filtered.
- GWN on a sphere: 1 inside, 0 outside to 1e-2; data-point query term counted.
- FWN vs direct GWN (3000-point ellipsoid, 2000 queries): β → ∞ equal to 1e-11; RMS by order at β = 2: 9.6e-3, 8.7e-3, **1.2e-3**; β = 4 better still. Single far cluster: error ratio between R = 4 and 8 is 0.154, 0.062, 0.026 for orders 0, 1, 2 (expected 2^{−(p+3)} = 0.125, 0.0625, 0.031). Triangles (torus mesh): max error < 1e-2, exact solid-angle sums are integers to 1e-9.
- Regularized WN: `S(1), S(2), S(3)` against tabulated erf to 1e-14, continuous at the series/CF switch, small-t asymptote, ε → 0 limit equals the GWN to 1e-10.
- SHC / SSPD / LSE gradients against central differences (constant λ with `All` and `Knn 32`, and `PaperHeuristic`) to 2e-6 relative; SSPD on a sphere ≈ SDF to 1e-2; shift invariance by a hand-written sum to 1e-13.
- Truth: cube mesh equals the analytic box SDF to 1e-14; icosphere (level 4) within the inscribed-ball bound `r − r_in`; torus mesh (128×64) max 4.5e-4 within the sagitta bound 4.5e-4; BVH closest point = brute force (distance bitwise); exact sign = FWN sign on 3000 queries; open / inward / bad index / repeated vertex / NaN errors; soups accepted with `require_closed: false`.
- Determinism: areas, SHC/SSPD/LSE values and gradients, GWN, FWN, truth SDF and winding bitwise equal at 1, 3, 8, 32 threads.
- Errors: non-finite queries, buffer size mismatches, negative / wrong-length areas, bad options, SHC `w = 0` → `DegenerateQuery` naming the query.

## Sanity run (WP14 preview)

`PAT_WP5_NAME=<input> PAT_WP5_THREADS=8 [PAT_WP5_SUBSAMPLE=512] target-pat5/release/deps/baselines-* sanity_grid_errors --ignored --nocapture`
(from the crate directory; 9 jobs in parallel, one log each in `runs_tmp/pat-wp5/sanity{,512}/`). Grid: the golden `grid_axis.npy`
(64³ over [−1, 1]³), mean absolute error against the truth. Truth: the ABC validation meshes (`ABC0_validation_0.npz` 1, 6;
`ABC17_validation_3.npz` 1, 3; `face_masks` applied, all closed, positive volume), the reference's normalized `block`/`knot` meshes,
`femur.obj`'s faces in the golden normalization (closed; the golden cloud is exactly its vertices), the analytic sphere/torus for
`synth_*`. Meshes extracted to `runs_tmp/pat-wp5/meshes/`. Baselines use Barill areas (k = 20) and the PAT neighbor rule (32-NN,
λ = 64/σ_x). "X-sign+nearest" = distance to the nearest cloud point signed by X (threshold ½); "LSE-GWN-signed" = the LSE soft-min
distance signed by the GWN. SPSR = Open3D depth 6, then `MeshTruth { require_closed: false }`. Our PAT-Reference grid equals the
authors' golden grid to ≤ 2.6e-7 on every input.

**Golden clouds as given** (1024–3897 points; synth 512):

| method | block | knot | femur | ABC0_1 | ABC0_6 | ABC17_1 | ABC17_3 | synth_sphere | synth_torus |
|---|---|---|---|---|---|---|---|---|---|
| PAT-Reference (32-NN) | 0.0172 | 0.0078 | 0.0018 | 0.0165 | 0.0367 | 0.0040 | 0.0135 | 0.0012 | 0.0026 |
| PAT-Geometric (32-NN) | 0.0172 | 0.0078 | 0.0018 | 0.0165 | 0.0367 | 0.0040 | 0.0135 | 0.0012 | 0.0026 |
| PAT-Reference (all) | 0.0177 | 0.0083 | 0.0032 | 0.0171 | 0.0369 | 0.0044 | 0.0142 | 0.0008 | 0.0028 |
| PAT-Geometric (all) | 0.0177 | 0.0083 | 0.0030 | 0.0171 | 0.0369 | 0.0044 | 0.0142 | 0.0008 | 0.0028 |
| SHC | 0.0145 | 0.0083 | 0.0109 | 0.0132 | 0.0389 | 0.0050 | 0.0128 | 0.0045 | 0.0083 |
| SSPD | 0.0589 | 0.0288 | 0.0349 | 0.0904 | 0.0887 | 0.0169 | 0.0859 | 0.0059 | 0.0190 |
| GWN-sign+nearest | 0.0126 | 0.0068 | 0.0001 | 0.0072 | 0.0265 | 0.0046 | 0.0069 | 0.0029 | 0.0048 |
| FWN-sign+nearest | 0.0126 | 0.0068 | 0.0001 | 0.0072 | 0.0264 | 0.0046 | 0.0069 | 0.0029 | 0.0048 |
| RWN(0.02)-sign+nearest | 0.0126 | 0.0068 | 0.0001 | 0.0072 | 0.0265 | 0.0046 | 0.0069 | 0.0029 | 0.0048 |
| LSE-GWN-signed | 0.0117 | 0.0060 | 0.0122 | 0.0068 | 0.0258 | 0.0055 | 0.0066 | 0.0036 | 0.0057 |
| SPSR | 0.0141 | 0.0061 | 0.0034 | 0.0130 | 0.0417 | 0.0046 | 0.0154 | 0.0013 | 0.0017 |

**512-point seeded subsets** (the §16.2 count; subset indices in `runs_tmp/pat-wp5/sub/*.idx`):

| method | block | knot | femur | ABC0_1 | ABC0_6 | ABC17_1 | ABC17_3 | mean (all 9) |
|---|---|---|---|---|---|---|---|---|
| PAT-Reference = PAT-Geometric | 0.0318 | 0.0118 | 0.0065 | 0.0253 | 0.0473 | 0.0127 | 0.0254 | 0.0183 |
| PAT (all) | 0.0317 | 0.0115 | 0.0070 | 0.0253 | 0.0470 | 0.0126 | 0.0254 | 0.0182 |
| SHC | 0.0384 | 0.0168 | 0.0154 | 0.0261 | 0.0673 | 0.0197 | 0.0255 | 0.0247 |
| SSPD | 0.0658 | 0.0364 | 0.0904 | 0.0888 | 0.1037 | 0.0157 | 0.0923 | 0.0576 |
| GWN/FWN/RWN-sign+nearest | 0.0313 | 0.0172 | 0.0050 | 0.0154 | 0.0373 | 0.0192 | 0.0190 | 0.0169 |
| LSE-GWN-signed | 0.0306 | 0.0160 | 0.0056 | 0.0144 | 0.0366 | 0.0185 | 0.0178 | 0.0165 |
| SPSR | 0.0259 | 0.0078 | 0.0079 | 0.0180 | 0.0590 | 0.0120 | 0.0230 | 0.0220 (7) |

Observations:
- **Reference vs Geometric:** identical on every input except full femur (the only one with saddle tori: 4), where Geometric is
  slightly better (32-NN 0.001775 → 0.001773; all 0.003217 → 0.003018). None of the 512-point subsets has a saddle torus. On
  `synth_torus` the network predicts no saddle even on the inner equator. A decision for WP14 needs clouds with saddles
  (holey-bunny has 33, but no mesh).
- PAT beats SHC and SSPD on most inputs, as in the paper. It does **not** beat "nearest-point distance signed by the GWN" or the
  GWN-signed LSE on these inputs, and it beats SPSR only on some (synth sphere, ABC0_6 at 512, femur). On femur the cloud is the mesh
  vertices, so the nearest-point distance is nearly exact (0.0001). The paper's "GWN" row may be defined differently; WP14 should
  state the definition used.
- SHC has `w = 0` exactly at grid points in the planes of axis-aligned CAD faces (ABC17_1: 143, ABC0_6: 167 of 262,144; 24 at
  512 points). The run sets φ = 0 there (`shc_zero_value: Some(0.0)`, counted); see proposal 4.
- FWN (β = 2) flips the GWN sign at 1–249 grid points per input (|Δw| up to 0.52 near the points). The MAE is unchanged to 4 digits.

## Decisions (implementation within the SPEC)

1. **Voronoi cell by half-plane clipping, not a planar Delaunay triangulation.** The cell of the center among its projected
   neighbors is `∩_j {y : ⟨y, q_j⟩ ≤ ‖q_j‖²/2}`, the dual cell of the center vertex in the Delaunay triangulation. Clipping a convex
   polygon by ≤ 20 half-planes needs no orientation predicates and no triangulation, and gives the same area.
2. **Unbounded cells:** every cell is intersected with the square of half-width `h = max_j ‖q_j‖` (the projected k-NN radius). Cells
   that touch it are counted as `clipped`. On a smooth closed surface none touch it (sphere test). Coincident projections share the cell.
3. `k = 20` neighbors **besides** the point (Barill: "each point and its k-nearest neighbors"). The reference's FWN path passes
   k = 21 including self (the same thing). The tangent plane is the point's normal plane (SPEC), not a best-fit plane (Barill).
4. **Own BVH in `truth.rs`.** It is not `ddg-repulsive-curves::Obstacle`. Reasons: no new dependency while another agent edits
   `Cargo.toml`; `Obstacle` rejects degenerate triangles, which CAD meshes and SPSR output can have; the deterministic tie rule on the
   face index. About 100 lines. The closest point follows Ericson §5.1.5, as `Obstacle` does.
5. FWN octree depth is capped at 48 (exact duplicates). Triangle leaves hold 8; point leaves hold 2 (Barill).

## SPEC change proposals

1. §9 point areas: replace "our own planar Delaunay (incremental with exact predicates …)" with decisions 1–2 above (same area,
   no predicates); add the clipping rule for unbounded cells and the `clipped`/`coincident`/`filtered` counters.
2. §13: the truth's closest point (Ericson 2004) and the triangle solid angle (van Oosterom–Strackee 1983) have no ids. I tagged
   the solid angle `pat-winding-dipole` (Barill Eq. 17: the integrated dipole) and left the closest point untagged. Proposal: new
   cited ids `pat-truth-closest-point` and `pat-solid-angle` in CITED.md, after which I retag.
3. §16.2 "GWN-sign + distance": define it. Ours: nearest-cloud-point distance signed by the GWN at ½. A GWN-signed LSE is also
   available. The table shows the definition matters.
4. WP3 resolution 2 (SHC `w = 0` → `DegenerateQuery`) makes batch SHC fail on CAD grids: exact zeros occur at in-plane grid points.
   Kept as the default. Added `KernelOptions::shc_zero_value: Option<f64>` (counted in `KernelReport::degenerate`) for error
   statistics. WP14 should fix the value it uses; 0 is the natural choice for in-plane points of a planar patch.
5. §9 FWN: Barill uses β = 2.3 for points (their RMS 1e-2 target) and 2 for triangles. The SPEC says 2 for both, which is the default here.

## Open issues

- `cargo clippy -p ddg-points-as-tori --all-targets -D warnings` currently fails on two `unused_mut` in `src/field.rs:521,536`. That
  is another agent's uncommitted edit, not WP5's.
- WP7 has its own `src/train/truth.rs`. If it computes mesh SDFs on the CPU, it could reuse `truth::MeshTruth`; the orchestrator
  should check.
- No performance numbers were tabulated (the reports carry `seconds`). The sanity jobs took 14–42 s each at 8 threads, including
  the brute-force GWN and RWN over 262k queries.
- SPSR meshes are closed except ABC0_6 (14 boundary edges at full density, 16 at 512 points: the surface leaves the
  reconstruction cube). Its volume is 1.47 against the true 0.32, so SPSR fails on this thin shape, and its GWN sign as a soup
  decides that row. Full femur has 3 non-manifold edges. Every SPSR volume is positive, so the paper's flip-if-negative check
  would not have triggered.
