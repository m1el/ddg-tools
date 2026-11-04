# WP2 notes (Lean proofs + citations)

Files: `lean/proofs/DdgProofs/PointsAsTori/*.lean`, `lean/proofs/DdgProofs/PointsAsTori.lean`,
the import line in `lean/proofs/DdgProofs.lean`, section "Points as Tori" of `lean/proofs/CITED.md`.
Namespace `DdgProofs.PAT`; vectors `RC.V3 = Fin 3 → ℝ`, `RC.vnorm`.
Build: `cd lean/proofs && lake build DdgProofs.PointsAsTori.<M> > runs_tmp/pat-wp2/<M>.log 2>&1`.

## Findings (to report; SPEC not edited)

1. **Saddle tori are mirrored (pat-torus-contact, Eq. 22 + center).** With the paper's center
   `c = o − sign(T)/|κmin| n*` and Eq. 22, the torus osculates `Q*` at `o` only when
   `κ₊κ₋ > 0`. For saddles (`κ₊κ₋ < 0`) the torus is the mirror image through the tangent plane:
   its curvatures at `o` w.r.t. `n*` are `−κmin, −κmax`, and its SDF has the wrong sign along
   `n*` (`g(o + εn*) = −ε`). Numerically (scratchpad `torus_check.py`): residual `g(Q*(h))` is
   `O(h²)` (1e-4 at h = 1e-2) with the paper rule, `O(h³)` with `c = o + n*/κmin`
   (the corrected center; the SDF sign of Eq. 22 is then right). Reference code
   `torus_distance.cpp:781-807` is the same as the paper. The pretrained network was trained with
   this rule, so inference must keep it; the SPEC claim "principal curvatures of the torus at `o`
   are `κmin, κmax`" holds only for `κ₊κ₋ > 0`.
2. **Zero eigenvector off the umbilic (pat-principal-pair).** `(v1 ± disc, v2)` is the zero
   vector for one sign whenever `v2 = 0` and `v1 ≠ 0` (e.g. `a10 = a01 = a11 = 0`, `a20 ≠ a02`).
   The reference clamps it to `(±1e-10, 0)` → the *other* eigenvector, so both curvatures become
   the same (`(a20, a02) = (−.5, −.2)` gives a sphere `κ = −1` instead of `(−1, −0.4)`). Near
   `v2 = 0` (|v2| ≲ 1e-8|v1|) the f64 cancellation in `v1 + disc` gives the same failure.
3. **Clamped exact zero curvature (pat-torus-contact with D9).** With `sgn(0) = 0`, an exact
   `κmin = 0` (cylinder, plane) gives `R = 1/ε`, `d = (0, 0)` at `o`: the torus does not pass
   through `o` (`g(o) = −sign·r`). Contact needs `κ₊κ₋ ≠ 0`.

## Progress

- [ ] Basic  - [ ] Features  - [ ] Fitting  - [ ] Torus  - [ ] Blend
- [ ] Baselines  - [ ] Loss  - [ ] Apps  - [ ] Pat2d  - [ ] CITED  - [ ] import line
