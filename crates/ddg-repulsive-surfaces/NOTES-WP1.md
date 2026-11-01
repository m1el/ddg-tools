# WP1 notes (crate skeleton, Surface, geometry, kernel, exact energy)

Build: `CARGO_TARGET_DIR=$PWD/target-wp1 cargo test -p ddg-repulsive-surfaces`.

## Done
- `Cargo.toml`, `README.md`, `LICENSES/repulsive-surfaces-MIT.txt` (copy of the reference's LICENSE).
- `src/lib.rs` declares every module of SPEC §15; stubs with module docs naming the WP that fills them.
- `error.rs` (`RsError`), `surface.rs` (`Surface`, `Component`, `Orientation`), `geometry.rs`,
  `kernel.rs` (`Exponents`, kernel + gradient, `ipow`), `energy/{mod,exact}.rs`, `io.rs`.

## Hooks for later packages (so they edit only their own files)
- `energy::energy` calls `energy::barnes_hut::energy_and_differential(pool, surface, positions,
  geom, exponents, opts) -> Result<(Ẽ, differential, pairs_exact, nodes_lumped), RsError>`
  (WP3 fills `barnes_hut.rs`; stub returns `RsError::Unsupported`).
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

## Remaining
(see below; updated as work proceeds)
