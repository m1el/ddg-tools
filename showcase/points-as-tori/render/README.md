# Sphere-traced renders (Fig. 12 analogues)

Made by `render.sh` (the `render` example of `ddg-points-as-tori`, CPU, 8 threads, 1024 x 1024,
reference weights, default evaluation options: 32-NN, per-query lambda). Shading is the authors' demo
shader (MIT, Copyright (c) 2026 Nicole Feng); see `crates/ddg-points-as-tori/NOTES-WP11.md`.

| image | scene |
|---|---|
| `tooth-level0.png`, `snail-level0.png`, `bunny-level0.png`, `analytic5-level0.png` | the zero level set |
| `tooth-offset-0.08.png`, `tooth-offset-0.16.png` | offset surfaces `phi = 0.08, 0.16` (`level=`) |
| `snail-offset-csg-0.05.png` | the offset as a CSG node (`offset=0.05`) |
| `boolean-subtract.png`, `boolean-union.png` | Booleans of two fitted fields (`op=`, `other=`) |
| `bunny-cut-plane.png` | cut plane `z = -0.1` coloured by phi with isobands |
| `tooth-demo-constants.png` | the demo's constants (safety 1, 40 steps, eps 1e-2) |

`safety-study.txt`: the tracer `safety` / `g_min` / step-budget study (SPEC 16.3) on three clouds.
Timings and ray statistics of every image: `runs_tmp/pat-wp11/showcase-render.log` (gitignored).
