# CCD investigation (SPEC §7.3, §16 "CCD default (D20)": EntangledDoubleTorus F = 42 with CCD vs 8.2 without)

Tool: `examples/ccd_investigation.rs` (`run`, `check`, `topo`, `klein`; see its header). Build:
`CARGO_TARGET_DIR=$PWD/target-ccdinv cargo build --release -p ddg-repulsive-surfaces --example ccd_investigation`.
`run` steps a scene like `run_scene --no-stop` (but `--no-ccd` turns CCD off, and CCD is forced on for
subcritical scenes) and checks every step end with `ccd::intersecting_pairs`. `check` needs the per-step
states (start, line-search trial, projection, remeshed end), which the public `Session` API does not
expose: the session was instrumented **locally and not committed** (patch at the end; env
`RS_CCDINV_DUMP=dir`, `RS_CCDINV_CAP=c`, `RS_CCDINV_NO_GUARD=1`). Logs, CSVs, dumps, OBJs:
`runs_tmp/ccdinv/` (not committed). All runs 2026-11-02, release, 8 threads each, 7 runs in parallel.

`check` per step: an exact swept test of the linear motion start → trial and trial → projection on the
same mesh (all vertex–triangle and edge–edge pairs without a common vertex, grid broad phase on swept
boxes; the coplanarity cubic's roots in [0, 1], a contact where the distance at the root is ≤ 1e-7 ℓ̄₀;
`side_change` = the cubic changes sign, i.e. an odd number of crossings); the library's
`first_contact` with gap 1e-4 and `close_factor 0` as a cross-check; a replay of the remesh pass from the
dumped projection (bit-identical to the session in every step of every run) with its largest vertex
displacement (vertex map), and the "sheet distance" (smallest distance between triangles more than 4
rings apart) before remeshing.

## Verdict: (A) — the run without CCD tunnels in step 6; (B) is negligible

* Steps 1–5 are identical with and without CCD (F 156.38 after step 5). In step 6 the start rule
  (`GrowPrevious`, doubling) proposes τ₀ = 3.72; without CCD the first trial is accepted, its farthest
  vertex moves **14.3 ℓ̄₀**. The swept test finds **2226 contacts (488 vertex–triangle, 1738 edge–edge),
  all with a side change**, at t ∈ [0.195, 0.588] of the trial: **244 vertices pass through 327 faces**,
  each vertex through two faces (into and out of a tube), around (0.07, −0.01, 1.80). First contact
  `EdgeEdge { a: 8633, b: 13175 }` at t = 0.1951 — the same pair and time the library CCD reports
  (`first_contact` gap 1e-4: τ_max = 0.1951; with the default gap the CCD run capped this step at
  ⅔·0.719 on the same `EdgeEdge { a: 8633, b: 13175 }`). The end of step 6 is intersection-free, so
  the step-end check of WP5 could not see it. Picture: the small ring of the double torus, which is
  linked through the large ring at the start, is pulled straight through the large ring's tube; after
  the step the two rings are unlinked.
* No other step of the no-CCD run (1–200) has a contact, in the line search or in the projection (≤ 1.2
  ℓ̄₀ moves). Remeshing: displacement ≤ 0.49 ℓ̄₀ per pass while distinct sheets are ≥ 2.0 ℓ̄₀ apart in
  every step, so no remeshing move can cross; every step end has 0 intersecting pairs.
* All CCD runs (default, gap 0.02, gap 0.01, cap 0.9, no guard, gap 0.01 + cap 0.9 + no guard): **0
  contacts in all 200 × 6 steps**, step ends intersection-free.
* Isotopy indicator: tree–cotree generators γ₁..γ₄ (genus 2 in every state), `V_ij = lk(γ_i, γ_j⁺)`
  (exact Gauss linking of polygons, push-off 0.05 ℓ̄ along the normal; entries integral to 1e-11) and
  `det(V − t Vᵀ)` (invariant under change of basis). It is `t²` for the initial mesh, both states
  around step 6 and all final states, and `det(V − Vᵀ) = 1`. This indicator does not separate the
  states (for a closed surface it is probably `t^g` always, since the kernels to the inside and the
  outside are complementary Lagrangians; not proved here). Topologically both end states are standard
  genus-2 surfaces: the start is the boundary of a neighbourhood of a handcuff graph with Hopf-linked
  loops, whose exterior (Hopf-link exterior T²×I minus the bar, an unknotting tunnel) is a handlebody,
  so by Waldhausen the surface is isotopic to the unlinked double torus. **So the no-CCD end state is
  (very likely) in the same isotopy class; its path is not an isotopy** (a regular homotopy through
  2226 crossings, which a barrier energy forbids and the paper's flow never takes). The fair
  comparison is the legitimate path, which has to slide one ring around the other through the bar.
* Shapes after 200 steps: CCD — two rounded Hopf-linked rings joined by a short bar, extents
  7.29 × 6.94 × 7.35, area 158.33 (constrained), volume 71.4, F = 41.8; no CCD — one rounded,
  unlinked genus-2 blob, 7.09 × 7.06 × 7.01, volume 109.2, F = 8.25.
* 700 steps (the scene's `iteration_limit`): see "Long runs" below.

## (B) measured: the CCD knobs barely matter here

EDT, 200 steps, all seven runs concurrently (same load):

| run | F(200) | steps capped by CCD | mean ms/step (1–140) | CCD ms/step | GMRES its/step (1–140) |
|---|---|---|---|---|---|
| CCD default (gap 0.05, cap ⅔, guard) | 41.80 | 3 (6, 8, 40) | 1954 | 17.1 | 34.3 |
| gap 0.02 | 37.27 | 8 | 1926 | 17.2 | 33.9 |
| gap 0.01 | 35.61 | 5 | 1930 | 17.4 | 33.8 |
| cap 0.9 | 40.67 | 4 | 1924 | 18.7 | 33.9 |
| no collapse guard | 41.80 (bit-identical: the guard never refused) | 3 | 1935 | 19.6 | 34.3 |
| gap 0.01 + cap 0.9 + no guard | 39.14 | 2 | 1926 | 16.5 | 34.1 |
| no CCD | 8.25 | — | 1456 | 0 | 21.9 |

The CCD caps only 2–8 of 200 steps; otherwise τ is set by the Armijo backtracking (≈ 0.94 halvings
per step in every run). The F differences between the CCD variants (35.6–41.8) are trajectory noise of
a few different early steps, not a systematic slow-down; they are all in the linked configuration.
CCD costs ~17 ms of ~1.9 s (0.9 %). The no-CCD run is faster per step because the unlinked shape
needs fewer GMRES iterations (22 vs 34), not because of the CCD. **Recommendation: keep the defaults**
(CCD on, gap 0.05 ℓ̄₀, cap ⅔, collapse guard). Smaller gaps or a larger cap buy nothing measurable here
and leave less margin; all variants were intersection-free along the whole path.

## Already-intersecting inputs (Klein bottle, `KleinBottle/scene.txt`, exponents (3.5, 7), subcritical)

200 steps, CCD forced on: F 793 → 631, τ ≈ 1e-4–1e-3, 149/200 steps capped; without CCD: F → 156
(as WP5). `klein` on the 200 no-CCD steps (library `first_contact`, default gap, on each step's motion):
the CCD would cap 193 steps (median τ_max 0.10 of the step) and in **191 of them the limiting pair
touches the current intersection set** (both primitives have a vertex of a face in
`intersecting_pairs` at t = 0, "ring 0"); 7 steps have no contact, 2 have a ring-1 contact. Of the
1011 exact crossings of the 200 motions, 919 are ring-0 pairs (the intersection curve sliding along the
mesh), 92 (9 steps, all after step 165 when the immersion simplifies) involve a pair one ring away, 15
(3 steps) two rings away, none farther.

Proposed rule (not implemented; `first_contact` would need an exclusion set): at t = 0 compute
`intersecting_pairs`, let `I` = the vertices of those faces, dilated by one ring; a vertex–triangle or
edge–edge pair is **unconstrained when both elements have a vertex in `I`**; all other pairs keep the
gap rule. On the Klein run this removes all 193 limiting contacts of the default CCD (191 ring 0, 2
ring 1); exact crossings outside the 1-ring band remain in 3 of 200 steps (15 pairs, all after step 185),
which the rule would cap — a limit on how fast the intersection curve may travel per step (2 rings: 0
steps). With it CCD can stay on for immersed inputs and
still forbids new intersections away from the existing curve. Until it exists, "CCD off for
subcritical exponents" (`scene_options`) is the right rule for subcritical scenes, which may intersect
by design; but the trigger should be "the input intersects" (`intersecting_pairs` non-empty at the
start) rather than the exponents, since a supercritical run on an immersed input stalls the same way
and a subcritical run on an embedded input loses its protection.

## Long runs (700 steps, the scene's `iteration_limit`; the two runs concurrently, no dumps)

The CCD run unlinks the rings legitimately (by sliding, no crossing at any step end; 5 capped steps in
700): F 41.8 (200) → 24.5 (230) → 14.9 (250) → 9.92 (273) → 8.43 (300) → **7.633 (700)**; without CCD
8.25 (200) → 7.630 (300) → **7.607 (700)**. Final shapes are the same unlinked genus-2 blob (extents
7.63 × 6.73 × 5.72 vs 7.44 × 7.03 × 5.89 up to rotation, volume 115.60 vs 115.48, `det(V − tVᵀ) = t²`).
So the "F = 42 vs 8.2 at 200 steps" gap is the ~70 steps it takes the legitimate flow to slide one ring
off the other, which the no-CCD run skipped by passing it through the tube in step 6. ms/step here
(run with 1 other job): 562 with CCD, 518 without (CCD 5 ms median).

## Local instrumentation used for the dumps (not committed)

```diff
@@ impl Session (step_inner, CCD cap)
             if rep.contact.is_some() {
-                tau0 = tau0.min(2.0 / 3.0 * rep.tau_max);
+                let cap: f64 = std::env::var("RS_CCDINV_CAP").ok().and_then(|s| s.parse().ok()).unwrap_or(2.0 / 3.0);
+                tau0 = tau0.min(cap * rep.tau_max);
             }
@@ (remeshing)
+            let no_guard = std::env::var("RS_CCDINV_NO_GUARD").is_ok();
             let mut guard = match &eff.line_search.ccd {
-                Some(o) => Some(CcdCollapseGuard::new(o, self.l0)?),
-                None => None,
+                Some(o) if !no_guard => Some(CcdCollapseGuard::new(o, self.l0)?),
+                _ => None,
             };
@@ (before "6. Commit.")
+        let ccdinv_dump = std::env::var("RS_CCDINV_DUMP").ok().map(|d| {
+            (d, self.surface.faces().to_vec(), self.positions.clone(), moved.clone(), proj.positions.clone())
+        });
@@ (after seconds.total)
+        if let Some((dir, f0s, p0, pm, pp)) = ccdinv_dump {
+            // little-endian: iteration u64, accepted u64, tau f64, |V0|, |F0|, |V1|, |F1| u64,
+            // F0 [u32;3]*, P0, trial, projection [f64;3]*|V0| each, F1, P1
+            ...
+            let _ = std::fs::write(format!("{dir}/step_{iteration:05}.bin"), b);
+        }
```

## Implementation of the proposed rule (branch of this worktree, 2026-11-02)

`ccd::IntersectionBand::new(pool, surface, positions)` runs `intersecting_pairs` (exact triangle test,
Morton box tree, parallel on the session pool) at the start of every step and marks the vertices of the
intersecting faces plus one ring (edge neighbours). `first_contact_with_band` (and `first_contact`, which
passes `None`) leaves a vertex-triangle or edge-edge pair unconstrained when both elements have a vertex
in the band (`CcdReport::exempt_pairs`, counted on the first window); mesh-obstacle pairs are never
exempt. `CcdCollapseGuard::set_band` skips a (new star triangle, triangle) pair when both have a band
vertex (a survivor counts as in the band if either end of its edge is; vertices created by the remesher
in the same pass are not). An empty band is `None`, so an embedded input runs exactly the old code path
(tests: empty band == `None` bitwise; `exempt_pairs == 0` at every step of an embedded flow).
`scene_options` no longer turns the CCD off for subcritical exponents; `RunSceneOptions::no_ccd` /
`run_scene --no-ccd` switch it off for comparison.

Klein bottle (`scene.txt`, 200 steps, `--no-stop`, 8 threads, remeshing on), objective `F`:

| run | F(200) | ms/step | GMRES its/step | CCD median share of a step |
|---|---|---|---|---|
| CCD on, with the rule | **25.29** | 97.8 | 34.0 | 2.2 % (includes the intersection test) |
| CCD off | 155.73 | 117.1 | 45.5 | - |
| CCD forced on, no rule (investigation) | 631 | - | - | - |

So the rule does not only remove the stall, it ends lower than the run without CCD (the unconstrained
run takes tiny Armijo steps in several steps: tau 1.6e-2 at step 200 against 2.9e-1). Sum of exempt
pairs over the 200 steps: 51,992.

Cost of the start-of-step test (`examples/ccd_band_cost.rs`, tori, best of 7, machine not otherwise
idle; `IntersectionBand::new` including the dilation; the CCD column is `first_contact_with_band` for
a push of 0.3 mean edge lengths):

| faces | input | 1 thread | 8 threads | CCD 1 / 8 threads |
|---|---|---|---|---|
| 10,000 | embedded | 19.4 ms | 7.3 ms | 46.8 / 16.4 ms |
| 25,600 | embedded | 36.6 | 12.3 | 111 / 22.6 |
| 50,000 | embedded | 69.7 | 16.5 | 216 / 80.7 |
| 20,000 | crossing (330 pairs) | 30.2 | 11.8 | 71 / 25.6 |
| 51,200 | crossing (522 pairs) | 78.9 | 28.7 | 215 / 61 |
| 100,000 | crossing (714 pairs) | 190 | 33.9 | 560 / 120 |

About 0.35 of the CCD itself, which is 1-2 % of a step, so ~0.5 % per step; there is no early exit to
take (proving that nothing intersects needs the whole tree traversal), the box test rejects almost all
face pairs and the common-vertex bookkeeping no longer allocates.

Tests: `tests/ccd.rs` (band empty for embedded meshes and bitwise equal to no band; the band frees a
sliding overlap but C-D contact of two disjoint spheres is found at the contact time of the pair alone;
random pushes of 3 of 4 bodies never leave a crossing outside the band at tau_max; threads 1/2/3/8
bitwise; guard refuses less with a band), `tests/flow.rs` (intersecting tori and embedded Hopf link,
threads 1 vs 3 bitwise, Klein bottle 200 steps `--ignored`, release, 6 min for both runs).

### Proposed SPEC text for 7.3 (append after "The collapse guard also refuses ...") and 16

> **Inputs that already intersect.** At the start of each step the faces that intersect (exact
> triangle-triangle test, BVH) are found; `I` is the set of their vertices widened by one ring. In
> the CCD a vertex-triangle or edge-edge pair whose two elements both have a vertex in `I` is
> unconstrained (counted in `CcdReport::exempt_pairs`); every other pair, including every pair between two
> regions that were disjoint at the start, keeps the gap rule. The collapse guard applies the same
> exemption to a (new star triangle, triangle) pair. For an embedded mesh `I` is empty and the CCD is
> bitwise unchanged. Mesh-obstacle pairs are never exempt. CCD is on for every exponent
> (`scene_options` no longer turns it off); an immersion is therefore free to slide its intersection curve
> (at most one ring per step) but cannot create a crossing elsewhere. Cost: the test is ~35 % of the CCD (7-16 ms on 10k-50k faces, 8 threads).
> Measured: Klein bottle, 200 steps, `F` = 25.3 (CCD on) against 155.7 (off) and 631 (forced on, no rule).
