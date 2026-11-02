# WP4 notes (dynamic remeshing)

Build: `CARGO_TARGET_DIR=$PWD/target-wp4 cargo test -p ddg-repulsive-surfaces --test remesh`
(debug ~30 s, release ~2 s; `--release -- --ignored timing --nocapture` for the 40k timing).

## Done
- `src/remesh/mesh.rs`: `EditMesh` (mutable mesh), `split` / `collapse` / `flip` with the SPEC §10
  single-operation numbering, refusals (`Refused { kind, message }`, `RefusalKind`,
  `RefusalCounts`), `CollapseGuard` / `CollapseQuery` (the CCD hook for WP8), `finish()` →
  `Remeshed { surface, positions, vertex_map, face_parents, pinned }` (revalidated).
- `src/remesh/geom.rs`: `circumcenter` (rs-circumcenter), `smoothing_displacement`
  (rs-smoothing-weights), `violates_delaunay` (rs-delaunay-test), `folds` (normal turn > π − 2 or
  degenerate), `FOLDOVER_LIMIT`.
- `src/remesh/transport.rs`: `VertexMap` (CSR new → old with weights): `entry`, `transport`
  (f64), `transport_vec3`, `transport_set` (threshold), `compose`, `identity`.
- `src/remesh/mod.rs`: `RemeshOptions` (+ `validate`), `RemeshOrder`, `EnabledOps`,
  `RemeshReport`, `remesh(surface, positions, pinned, l0, options, guard)`.
- `tests/remesh.rs`: the 16 oracle `remesh-ops` fixtures (bitwise positions, faces, parents,
  vmaps, pins); every refusal kind; 24 stress runs (6 meshes × 4 seeds, threads) of 1500 random
  ops with `finish()` (full `Surface` + position revalidation) after **every** op; 4 full passes ×
  2 anisotropic scalings on every mesh; stretched icosphere quality; pins; boundary; flat-patch
  Delaunay; determinism (both orders); transport of linear data; 40k timing (ignored).

## Measurements (release, 1 thread)
- Stretched icosphere(3) scaled (4, 1, 0.5), L₀ = 1.5 × the unscaled mean edge: edges outside
  [L₀/2, 3L₀/2] 43.5 % → 0.00 % after 8 passes (one pass already splits 784, collapses 426);
  min angle 9.19° → 25.3° (1 pass) → 32.3° (8 passes); 1st percentile 11.0° → 35.2°.
- Stress: 24 runs × 1500 random ops (≈ 10k splits, 8k collapses, 7.5k flips, 10k refusals,
  reordered flips on Möbius/Klein), valid after every op.
- One pass on torus 200×100 (40 000 faces) scaled (1.3, 1, 0.8): 0.20 s (11k splits, 3.5k
  collapses, 11.7k flips → 54 882 faces); steady state on 48 000 faces: 0.13 s.
  (Profile: smoothing ≈ ½ — fans and the angle-weighted normals; flips ≈ ¼.)

## Decisions
- **Representation.** Not a pointer-based half-edge structure: face slots + per-vertex face lists;
  a half-edge is (face, side), its twin found through the vertex lists. This is orientation-free
  (non-orientable surfaces need nothing special except the flip's local reordering and the
  consistently oriented fan for the vertex normal). Dead slots are dropped in `finish()`, keeping
  relative order, which gives the oracle's numbering for single ops.
- **Foldover** = the normal of the *same* face (same vertex order) turns by more than π − 2, or the
  face becomes degenerate (SPEC/oracle). The reference instead measures the dihedral angle between
  each new face and the face across its link edge (`ref: remeshing.cpp:108–118, 173–185`).
- **DegenerateFace** (new refusal kind): split, collapse and flip always refuse a zero-area result
  (found by the stress test: splits of nearly degenerate faces and unchecked collapses/flips), so
  `finish()` always revalidates. With the foldover check on, a degenerate result is `Foldover`
  (as the oracle).
- **Oracle `InconsistentOrientation`**: the Rust flip reorders and reports
  `FlipOutcome::reordered`; the fixture test asserts that instead.
- **Smoothing** directions are Jacobi (all from the positions at the start of the smoothing);
  they are applied in vertex order, each guarded against the face normals *at the start of the
  smoothing* with the current positions of the neighbours, so no face turns by more than π − 2 in
  a round. Vertices with no fan / no normal are skipped (`smoothing_skipped`).
- **Collapse order** in the pass: shortest first (heap, ties by vertex ids), survivors' new short
  edges queued again, stale entries skipped. Splits: input edges in lexicographic order.
- **Flips**: FIFO of all edges (lexicographic), the 4 quad edges re-queued after a flip; the cap
  `flip_cap_per_vertex · |V|` per round; `flip_cap_hits` counts rounds that hit it.
- **Delaunay test** in the cotan form `d_c s_d + d_d s_c < 0` (no `atan2`), the foldover test as
  `before·after < cos(π − 2)|before||after|`; both are the tagged equivalences, and replay the
  oracle fixtures.
- **Options**: `split > 2 · collapse` required (halves of a split edge must not be collapsible);
  `rho ∈ (0, 1]`; `flip_cap_per_vertex ≥ 1`; `l0` finite > 0.
- **Pins**: `Remeshed::pinned[i]` is the new id of input pin `i` (same order). A collapse into a
  pinned vertex keeps the pin's position; its vertex-map entry has weight 1 on the pin's old id
  (the other end listed with weight 0).
- **Errors**: a failed revalidation returns the `Surface::new` / `validate_positions` error as is
  (there is no `RsError::Internal`); it has not occurred in any test.

## Remaining / for WP5 and WP8
- WP8: implement `CollapseGuard` (CCD gap test) and pass it to `remesh`.
- WP5: see "How to call" below.

## How to call (WP5)
```rust
let (out, report) = remesh(&surface, &positions, &pins, l0, &opts, None)?;
// out.surface, out.positions: the new mesh; out.pinned: new ids of the pins (same order)
// per-vertex data: out.vertex_map.transport(&f64s) / transport_vec3(&targets) / transport_set(&set, 0.5)
// out.face_parents[f]: old face of new face f (Blender material index)
```
`l0` = mean edge length of the run's **initial** mesh (`geometry::mean_edge_length`), fixed.
