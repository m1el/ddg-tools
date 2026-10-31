# ddg-tools

Rust implementations of discrete differential geometry papers, with Python bindings and a Blender
add-on built on them. Every formula the libraries use is derived in Lean (proved with Mathlib where
it allows, cited otherwise), and a Lean implementation serves as a differential-testing oracle.

## Papers implemented

- Yousuf Soliman, Dejan Slepčev, Keenan Crane. *Optimal Cone Singularities for Conformal Flattening.*
  ACM Transactions on Graphics 37(4), 2018. https://doi.org/10.1145/3197517.3201367
  (`crates/ddg-optimal-cones`)
- Rohan Sawhney, Keenan Crane. *Boundary First Flattening.*
  ACM Transactions on Graphics 37(1), 2017. https://doi.org/10.1145/3132705
  (`crates/ddg-bff`: flattening with cones, and conformal mapping of a genus-0 mesh to the sphere)
- Nicholas Sharp, Keenan Crane. *Variational Surface Cutting.*
  ACM Transactions on Graphics 37(4), 2018. https://doi.org/10.1145/3197517.3201356
  (`crates/ddg-variational-cuts`). A clean-room implementation from the text of the paper.
- Nicholas Sharp, Yousuf Soliman, Keenan Crane. *Navigating Intrinsic Triangulations.*
  ACM Transactions on Graphics 38(4), 2019. https://doi.org/10.1145/3306346.3322979
  (`crates/ddg-intrinsic`, and the applications of its Sec. 5 in `crates/ddg-intrinsic-apps`)
- Oded Stein, Eitan Grinspun, Keenan Crane. *Developability of Triangle Meshes.*
  ACM Transactions on Graphics 37(4), 2018. https://doi.org/10.1145/3197517.3201303
  (`crates/ddg-developability`). A clean-room implementation from the text of the paper.
- Chris Yu, Henrik Schumacher, Keenan Crane. *Repulsive Curves.*
  ACM Transactions on Graphics 40(2), 2021. https://doi.org/10.1145/3439429
  (`crates/ddg-repulsive-curves`). Parts are ported from the authors' MIT code, with its licence
  notice.

Each crate's README credits the further papers it uses, and its `SPEC.md` (where there is one) is
the discrete method as implemented, listing every difference from the paper as printed.

## Layout

| Path | What |
|---|---|
| `crates/ddg-core` | mesh topology and validation, geometry, cotan operators, shared numerics |
| `crates/ddg-optimal-cones` | the optimal-cones solver |
| `crates/ddg-bff` | boundary first flattening, cutting to a disk, mapping to the sphere |
| `crates/ddg-variational-cuts` | the cut flow (two regions and networks of regions) |
| `crates/ddg-intrinsic` | signposts, flips, Delaunay flipping and refinement, optimal Delaunay, overlay, point queries |
| `crates/ddg-intrinsic-apps` | Steiner trees, FEM conditioning, heat-method distance, adaptive refinement, tangent fields |
| `crates/ddg-developability` | the six developability energies, L-BFGS, subdivision, remeshing, cut and flatten |
| `crates/ddg-repulsive-curves` | the tangent-point energy, fractional Sobolev descent, constraints and penalties, collision-guarded steps, L-BFGS/NCG, Barnes-Hut, hierarchical matrices and multigrid for curve networks |
| `crates/ddg-uv` | end-to-end UV pipeline: polygons to charts, cones, cuts, flattening and per-corner UVs |
| `crates/ddg-py` | PyO3 bindings (Python module `ddg_tools`) |
| `blender/` | the Blender add-on: Conformal Unwrap (Optimal Cones); Cut: Start, Flow, Apply; Conformal Map to Sphere; Make Developable, Show Developability, Developable Seams + UV; Repulsive Flow (Live, Steps) |
| `lean/proofs`, `lean/oracle` | the derivations (`CITED.md` lists what is cited rather than proved) and the oracle |
| `showcase/` | stills and animations of every capability, rendered from the libraries' own results |
| `scripts/` | builds of the wheel and the Blender extension, test runners, `check-formulas.sh` |

See `CLAUDE.md` for the conventions shared by all crates and the mathematical contract.
`ROADMAP.md` lists what is scheduled next and the known gaps.

## License

MIT. Ported code keeps its original license notice in the crate that contains it.
The datasets used, with sources, citations and terms, are listed in [`DATASETS.md`](DATASETS.md).
