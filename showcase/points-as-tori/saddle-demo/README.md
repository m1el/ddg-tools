# Points as Tori: saddle demonstration

Figures for `crates/ddg-points-as-tori/SADDLE-DEMO.md`. It covers saddle-dominated analytic shapes:
a fat torus (R .54, r .36) sampled uniformly, mostly on its inner half, or only on its inner half;
and a hyperboloid-waisted spool. Each shape is compared under four methods: the authors' released network
(N), the exact fundamental forms through the paper's torus fit (E-ref), the same with the one-line
centre fix / our `SignRule::Geometric` (E-fix / E-geo), and a GWN-sign + nearest-distance baseline.

| file | made with | what |
|---|---|---|
| `render_<case>.png` | **our** renderer (`ddg-points-as-tori` CPU sphere tracer) on our crate's fields, which equal the authors' to < 1e-6 | columns truth / N / E-ref / E-geo; rows: surface, cutaway (camera-side half removed), cut plane y = 0 coloured by φ (isobands every 0.05, blue inside) |
| `slice_<case>.png` | matplotlib from **the authors' program's** arrays | φ on the plane y = 0 through the axis (solid: the method's zero set, dashed: truth) and the error map |
| `tori_at_saddles.png` | matplotlib from **the authors' program's** tori | individual tori at three saddle points, cut by the two principal planes: truth, network torus, paper-rule torus, fixed torus |
| `band_mae_by_region.png` | matplotlib from **the authors' program's** grids | narrow-band MAE split by the sign of K, all eight cases |

Cases: `torus-uniform-{512,2048}`, `torus-inner80-{512,2048}` (80 % of the points on the inner,
saddle half), `torus-inneronly-{512,2048}` (open patch), `spool-{512,2048}`. Renders are 384² per panel,
montages ≤ 1024 wide. Inputs, drivers and metrics: `crates/ddg-points-as-tori/tools/saddle_demo/`
and `examples/saddle_demo.rs`.
