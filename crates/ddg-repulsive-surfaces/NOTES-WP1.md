# WP1 notes (crate skeleton, Surface, geometry, kernel, exact energy)

Build: `CARGO_TARGET_DIR=$PWD/target-wp1 cargo test -p ddg-repulsive-surfaces`.

## Done
- `Cargo.toml`, `README.md`, `LICENSES/repulsive-surfaces-MIT.txt` (copy of the reference's LICENSE).
- `src/lib.rs` declares every module of SPEC §15; stubs with module docs naming the WP that fills them.
- `error.rs` (`RsError`), `surface.rs` (`Surface`, `Component`, `Orientation`), `geometry.rs`,
  `kernel.rs` (`Exponents`, kernel + gradient, `ipow`), `energy/{mod,exact}.rs`, `io.rs`.

## Hooks for later packages (so they edit only their own files)
- `energy::energy` calls `energy::barnes_hut::energy_and_differential(pool, surface, positions,
  geom, exponents, opts) -> Result<EnergyEvaluation, RsError>` (`EnergyEvaluation { energy,
  differential, pairs_exact, nodes_lumped }` in `energy/mod.rs`; WP3 fills `barnes_hut.rs`; the
  stub returns `RsError::Unsupported`).
- `energy::energy` calls `penalties::evaluate(pool, surface, positions, penalties, &mut
  differential) -> Result<Vec<f64>, RsError>` (weighted values; adds weighted gradients). WP7
  replaces the empty `Penalty` enum and the hook in `penalties/mod.rs`.
- `Energy { TangentPoint, Coulomb, Willmore }` is declared in `energy/mod.rs`; `energy()` evaluates
  the tangent-point energy only. Dispatch on `Energy` belongs to the session (WP5) / WP7.
- `BhOptions { theta, leaf_size }` (default 0.5, 8) and `validate()` are in `energy/mod.rs`.
- Error variants for later WPs are declared in `error.rs` (fields guessed; extend if needed —
  that is an edit of a WP1 file, coordinate via the orchestrator).
- `geometry::{face_pull_back, pull_back}` turn per-face partials `(∂/∂a, ∂/∂X, ∂/∂N)` into the
  vertex differential (WP3 can reuse for the BH gradient).

## Decisions
- `Surface::orientation(c)` returns `&Orientation` (SPEC §9.1 says `Orientation` by value; the
  `Orientable` variant carries a Vec).
- `Orientable { flips }`: the component's lowest face keeps its orientation.
- `Component::inconsistent_edge()` gives the witness edge for `VolumeNeedsOrientation`.
- `EnergyReport.differential` is the differential of `energy + Σ penalties`; `energy` excludes
  penalties (SPEC underspecified).
- `RsError::NonFiniteEnergy { faces: [u32; 2] }`: first non-finite pair in row-major order, else
  the pair with the largest term (overflow of the sum).
- `RsError::PolygonFace { face, corners }` added for OBJ polygons; `Unsupported { what }` for stubs.
- `read_obj` does not validate positions (callers do).

- `energy::exact::{energy_value, face_partials, energy_and_differential, make_pool}` are public
  (value-only for line searches; per-face partials for reuse).
- `kernel::ipow` (binary powering) is the integer path: `α` integer ≤ 64, `β` integer ≤ 128
  (`r^{-β}` from `r²`, one `sqrt` for odd `β`); others `powf`.
- `[formula: rs-kernel-gradient]` is tagged in `kernel.rs` (2 places). SPEC §12 allows WP10 to
  drop the Lean theorem ("FD + oracle, no tag"); if so, remove those two tags.

## Tests (all pass; clippy --all-targets clean; rustdoc -D warnings clean)
- `tests/surface.rs`: every topology/position error class with messages; icosphere vs `TriMesh`
  (edges, fans); Orientable flip sets; Möbius band (NonOrientable, 1 loop of 24); the reference's
  `klein_bottle_blocky.obj` (read from its path; NonOrientable, χ = 0, closed); 3 components with
  loops vs `TriMesh` (same cyclic order); `W`, `M` vs `ddg_core::operators` + flip invariance;
  face derivative FDs; OBJ round trip, polygon and parse errors.
- `tests/energy.rs`: exponent validation and identities; `ipow` vs `powf`; kernel FD + ν = 0
  continuity; energy FD on 5 meshes × 3 exponent pairs (rel. errors 1e-10..7e-9, run in
  parallel); rigid/scaling invariance and the three gradient identities; face flips/rotations;
  two hand-computed two-triangle energies (one bit-exact); naive double loop; bitwise equality at
  1, 3, 8 threads; NonFiniteEnergy / InvalidOption errors.

## Remaining
Nothing of WP1. Resume: run the tests above; see the hooks list for what later WPs plug in.
