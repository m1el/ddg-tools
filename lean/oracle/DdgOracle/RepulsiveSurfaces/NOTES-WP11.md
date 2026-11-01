# WP11 notes: the `rs` differential-testing oracle

Contract: `crates/ddg-repulsive-surfaces/SPEC.md` §13.5. Written from the SPEC only. The protocol
(fields per task) is in `lean/oracle/README.md`, section "Repulsive surfaces".

## Done

* Tasks `energy`, `df`, `metric`, `constraints`, `saddle`, `precond`, `step`, `remesh-ops`
  (`Energy`, `Metric`, `Constraints`, `Saddle`, `Step`, `Remesh`, `Format`).
* `rs-gen` (`Gen.lean`: tetra, octa, icosa, torus grids, patch with boundary, Moebius band, two
  components; perturbed), `rs-selftest [N]` (`Selftest.lean`; 17k checks at N = 120, 0 failures),
  `rs-fixtures write|check DIR` (`Fixtures.lean`, 70 cases in `lean/oracle/fixtures/rs/`,
  `MANIFEST.tsv` = `name task ok|error description detail`).
* Main dispatch (`rs`, `rs-gen`, `rs-selftest`, `rs-fixtures`), README section.

## Remaining (not WP11)

* `crates/ddg-repulsive-surfaces/tests/oracle.rs` (replay + `#[ignore]` fuzz), after the crate exists.
* Not modelled: penalties, Barnes-Hut, hierarchical products, GMRES, CCD, `BarycenterMode::Global`,
  `GrowPrevious`/`Reference` Armijo, the other metrics (§4.5), obstacles, the full remeshing pass.

## How to resume

```
cd lean/oracle
lake build > ../../runs_tmp/x.log 2>&1
.lake/build/bin/ddg-oracle rs-selftest 120
.lake/build/bin/ddg-oracle rs-fixtures check fixtures/rs      # tolerance 1e-9 relative to max(|x|,1)
.lake/build/bin/ddg-oracle rs-fixtures write fixtures/rs      # after an intended change; review the diff
.lake/build/bin/ddg-oracle rs-gen SEED COUNT [FIRST] | .lake/build/bin/ddg-oracle rs
```

## Readings and choices where the SPEC leaves room

1. **Newton correction sign (SPEC issue).** §6.1 defines `S = −C_N A′⁻¹ C_Nᵀ` and §6.2 says
   `h = A′⁻¹ C_Nᵀ S⁻¹ Φ` with `C_N h = Φ`, subtracted. With that `S`, `C_N h = −Φ`. The oracle uses
   `h = −Z S⁻¹ Φ` (`Z = A′⁻¹ C_Nᵀ`, position block), which gives `C_N h = Φ`, so `f ← f − h` removes
   the violation. The direction `x = y + Z S⁻¹ C_N y` (`y = A′⁻¹ g`) is correct as printed
   (`C_N x = 0`). Proposed SPEC text for §6.2: "`h = −A′⁻¹ C_Nᵀ S⁻¹ Φ(f)` (with `S` as in §6.1, so that
   `C_N h = Φ`) is subtracted".
2. **Newton with the old Jacobian can diverge (SPEC remark).** Area and volume gradients are nearly
   parallel on near-spheres (`∇A = (2/r)∇V` on a sphere), so `S` is nearly singular for both rows
   together; an accepted step that moves vertices by a mean edge length with targets 3 % off made the
   frozen-Jacobian iteration blow up to 1e80 in the selftest. §6.2 should say what the Rust code does
   when `|Φ|` grows (stop and report). The oracle simply iterates `newton N` times and prints
   `phi-before` / `phi`; it never reports a failure.
3. **Barycenter target with `component: None`.** §5.1 allows `target: Some(x)` with all components.
   The oracle applies the same target to every component. Suggested SPEC text: "a target with
   `component: None` applies to every component; use one constraint per component for different targets".
4. **Component numbering** (needed by `Barycenter`/`Area`/`Volume` with a component index): by smallest
   face index. Suggested SPEC text in §2.1.
5. **Redundant rows (D10).** The oracle never drops rows: a rank-deficient `A′` (or a singular `S`)
   is `SingularSaddle` (equilibrated pivot ratio below `1e-13`). The Rust crate drops redundant rows
   and reports; fixtures and generated cases avoid redundancy (the error fixture `rs-error-singular`
   pins the same vertex twice, which the Rust crate may legitimately repair; the replay test should
   skip it or expect either outcome).
6. **`precond`** takes positional rows only (nonlinear rows: `Unsupported`), prints `L′` with the
   rows in constraint order. `L_{2−s}` is `L_τ` of §4.2 with `τ = 2 − s`.
7. **`step`.** Gradient `g = dℰ̂` only (no penalties). Start rules: `abs T`, `maxdisp L`
   (`τ₀ = L ℓ̄₀ / max_v |x_v|`), `reference` (`2/‖x‖_F`). Halve until Armijo or `halvings > maxhalv`
   or `τ < 1e-12 ℓ̄₀/‖x‖_∞` (stuck: positions unchanged, projection still runs). Projection: positional
   exactly (barycenter translation with current dual areas, pins reset), then Newton rounds
   (`|Φ_r| ≤ 1e-8 |target|`, or `1e-8` for target 0, else `newton fixed N`: exactly N rounds), each
   followed by the positional projection. The Newton rounds reuse `Z`, `S` of the step (old Jacobian).
8. **`remesh-ops`** numbering (fixed in `Remesh.lean`): split appends the midpoint as vertex `NV` and
   the new face of each cut face is appended in increasing parent order; collapse keeps `min(a,b)`,
   removes `max(a,b)` and shifts higher indices down; position: pinned endpoint, else the single boundary
   endpoint, else the midpoint (both boundary: midpoint; both pinned: `BothPinned`). Flip: lower-index
   face `F1` traverses `p → q`, `F1 := (p, d, c)`, `F2 := (d, q, c)`. Always checked: edge exists, link
   condition (common neighbours = incident faces; an interior edge between two boundary vertices
   fails), no repeated vertex triple afterwards, the result is a valid manifold. Optional: `legality
   delaunay` (opposite angles sum `> π`), `legality foldover` (a surviving/new face turning its normal
   by more than `π − 2` or becoming degenerate). Flips need consistently traversed faces
   (`InconsistentOrientation` otherwise) and no existing edge `c d` (`FlipNonManifold`).
9. **`metric`'s `B⁰`** is built with the ordered kernel `K₂(S,T) ρ^{−2s}` directly (the symmetric form
   is checked against it in the selftest). `W`, `M` are orientation free.
10. Kernel `|ν|^{α−2}ν` is evaluated as `sign(ν) |ν|^{α−1}`; `ν = 0` gives 0.
