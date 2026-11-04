# WP14 notes: evaluation and black-box comparison

Results: `EVALUATION.md` (accuracy, continuity, sign rule, orientation, density, recommendations) and
`COMPARISON.md` (ours vs the authors' program). This file holds what the orchestrator needs to act on:
SPEC change proposals, open items, and how to resume.

## Files

| file | content |
|---|---|
| `examples/eval.rs` | `curv` (predicted vs true principal curvatures), `prepare` (datasets), `case` (64³ errors of every method, black-box comparison, slices), `continuity` (jump statistic, error, eikonality, MC components per rule at 64/128/256, mechanism at the top edges), `orient` (Fig. 27 at 10 and 30 iterations), `swaps` (size of the real kNN-swap discontinuities), `profile` (width of the steep transitions), `gpu` (TF32 vs Split3, needs `--features cuda` and a free GPU), `selfint` (the Thingi10K self-intersection test on single meshes) |
| `examples/density.rs` | Fig. 7 / Fig. 29 study |
| `examples/eval_util/util.rs` | npy writer, grids, metrics, jump statistic, JSON writer, `TaperField` (tapered kNN, ours, example code only), self-intersection test |
| `scripts/pat_eval.sh` | the whole pipeline by stage |
| `scripts/pat_eval_data.py` | ABC extraction from the training data and the train/test overlap check |
| `scripts/pat_eval_ref.py` | the authors' program on our clouds (their venv, CPU JAX) |
| `scripts/pat_eval_report.py` | tables and figures from the per-case JSON |

Work directory (gitignored): `runs_tmp/pat-wp14/` (`cases/<dataset>/<case>/` with `points.npy`,
`normals.npy`, `meta.txt`, truth mesh, `spsr/`, `ref/`, `result.json`, `slices.npy`,
`continuity.json`, `orient.json`; `report/tables.md` and `report/fig/`; logs in `logs/`).

## SPEC change proposals

1. **§16.2 "100 ABC shapes not seen in training" cannot be met with the shipped data.** The test splits
   of chunks 2/17/55/58 are not shipped. The only ABC shapes outside the training files are the 20 of
   `ABC0_validation_0.npz` and `ABC17_validation_3.npz`, and 2 of those are exact duplicates of training
   shapes (`abc0val_3` = `nice/ABC2_train_mode0_4.npz:205`, `abc17val_0` = `nice/ABC17_train_mode0_1.npz:131`;
   sorted vertex coordinates equal to 1e-6; `runs_tmp/pat-wp14/src/overlap.json`). WP14 uses the 20 plus
   80 seeded picks from the training files and reports the three groups (18 unseen, 2 duplicates, 80 seen)
   separately. Proposal: §16.2 says so. The 100 blobs are new (our generator, WP8), so unseen by
   construction.
2. **§16.2 "GWN-sign + distance"** (WP5 resolution 3) is reported in every table with two more criteria:
   narrow-band MAE (`|φ_true| < 0.05`, about 1.6 grid cells of 64³) and eikonality (mean `| ‖∇_h φ‖ − 1 |`
   by central differences on the grid, the same for every method, so methods without a gradient are
   comparable). Proposal: add these definitions to §16.2, plus "wrong sign on more than 1 % of the grid
   points with `|φ_true| > 0.05`" as the definition of a catastrophic failure (the paper does not define it).
3. **WP6 resolution 1 (the released blend is discontinuous): the mechanism is different.** See
   EVALUATION.md "Continuity". Proposal for §7.4/§8.6: "kNN and radius neighborhoods make φ discontinuous
   where the set changes, but the jumps are small (measured ≤ 1.5e-3, median ≤ 3e-5, at sampled set
   changes). The large `|Δφ|/h` that grows like `1/h` on grids down to `h ≈ 1/λ_x = σ_x/C` is a continuous
   transition of width ≈ 1–2/λ_x where the soft-argmin over data points switches between points whose tori
   disagree (`Δg` up to 0.5). Its slope is ≈ `λ_x Δg` (measured up to 90). The `All` rule has `σ_x` ≈ half
   the cloud diameter, hence smaller `λ_x` and bounded slopes." The tapered kNN of the WP6 proposal does
   not change the measures (it removes only the tiny swap jumps) and is not recommended.
4. **New finding for §7.1/§10 (the released weights): the tori are near-spheres.** At every input point
   of the 70 App. E shapes, the predicted |κmin|/|κmax| has a median of 0.89–0.95 (true value 0 on
   cylinders and extrusions). On saddle surfaces (35–69 % saddle area) the network predicts no saddle
   at all. The authors' golden tori are 97.6–100 % spindle tori (R < |r|), 81–100 % with R < |r|/4.
   So `SignRule` matters only on the rare predicted saddles (51 of 258,478 points), where Geometric
   is better on every cloud. Proposal: §18 records that the released network is effectively
   "points as spheres", that the `SignRule` default can follow geometry (Geometric) without
   measurable risk, and that WP10's comparison should report the anisotropy and saddle statistics
   of the retrained model (`eval curv <case>`).
5. **§16.3 "Tracer safety" was not measured (WP11's tracer; a timed study for WP16). ε was measured by
   WP9a (flat); not repeated.

## Open items

- GPU `GpuMatmul::Tf32` vs `Split3` effect on φ: the training run held the GPU. The `gpu` subcommand is
  ready: `cargo build --release -p ddg-points-as-tori --example eval --features cuda`, then
  `eval gpu <case-dir>... --threads 8` (writes `<case>/gpu.json`). For WP16.
- Retrained rows: `RETRAINED="--retrained retrained=data/points-as-tori/retrained-<name>.safetensors" scripts/pat_eval.sh case report`
  (append `:geometric` to the path for weights trained with `SignRule::Geometric`).
- Fig. 29 stops at 7.3M points (29M needs ~4 h of CPU network time on a shared machine).
- `DATASETS.md` (repo root, not WP14's file) should credit the Thingi10K models named in EVALUATION.md;
  EVALUATION.md has the credit table.
