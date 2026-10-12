# ddg-core

Foundation crate of the `ddg-tools` workspace: what every discrete differential geometry
algorithm here needs before it can start.

| module        | contents |
|---------------|----------|
| `topology`    | `TriMesh`: validated, oriented, manifold triangle mesh with halfedge connectivity, unique edges, boundary loops, components, Euler characteristic and genus; explicit repair helpers `orient_faces` and `compact_vertices` |
| `triangulate` | `triangulate_polygons`: splits polygons exactly like Blender and maps every triangle corner back to its polygon corner (cargo feature `blender-tessellation`, on by default) |
| `geometry`    | `Geometry`: corner angles, areas and cotangents from edge lengths (from positions, or intrinsic) |
| `operators`   | cotan Laplacian, lumped mass (raw and area normalized), angle defect, boundary geodesic curvature, Gauss–Bonnet residual |
| `sparse`      | `SparseMatrix` (CSC) and direct solvers `LLᵀ`, `LDLᵀ`, `LBLᵀ`, `LU` with reusable symbolic analysis, explicit orderings and an explicit thread pool |
| `io`          | Wavefront OBJ reader and writer (polygons, per-corner UVs), cone files of the optimal-cones dataset |
| `generators`  | icosphere, torus, flat and bumpy grid disks, annulus, open cylinder, seeded perturbation |

Conventions (indices, signs, normalization) are documented in the module docs and follow the
repository's `CLAUDE.md`. In short: halfedge `3f + k` runs from `faces[f][k]` to
`faces[f][(k+1) % 3]`, and the same id names the corner of face `f` at `faces[f][k]`; edges are
`[min, max]` pairs in lexicographic order; the Laplacian is positive semidefinite; the angle
defect lives on interior vertices and the geodesic curvature on boundary vertices. Invalid input
is reported by error enums that say what is wrong and where; nothing panics, and nothing is
repaired, normalized or clamped unless the caller asks for it.

A mesh with a vertex that no face uses is rejected (it would give the Laplacian an empty row);
`compact_vertices` removes such vertices and returns the index maps.

## References

- Yousuf Soliman, Dejan Slepčev, Keenan Crane. *Optimal Cone Singularities for Conformal
  Flattening.* ACM Transactions on Graphics 37(4), 2018.
  <https://doi.org/10.1145/3197517.3201367> — the discrete operators of Sec. 5.1.
- R. H. MacNeal. *The Solution of Partial Differential Equations by Means of Electrical
  Networks.* PhD thesis, Caltech, 1949 — the cotan operator.
- W. Kahan. *Miscalculating Area and Angles of a Needle-like Triangle.* 2014 — the numerically
  stable area and angle formulas used in `geometry`.
- Blender, <https://www.blender.org> — `triangulate` is a port of its face tessellation
  (see License). Its 2D polygon fill is in turn based on libGDX's `EarClippingTriangulator`.
- faer, <https://codeberg.org/sarah-quinones/faer> — the sparse factorizations behind `sparse`.

## Formula tags

Functions that implement a mathematical formula carry a `[formula: <id>]` line in their docs.

| id | function | statement |
|----|----------|-----------|
| `euler-characteristic` | `TriMesh::euler_characteristic` | `χ = V − E + F` |
| `genus-from-euler` | `TriMesh::genus` | `g = (2C − b − χ) / 2`: total genus of an orientable surface with `C` components and `b` boundary loops |
| `edge-length-euclidean` | `geometry::edge_length` | `ℓ = √(Δx² + Δy² + Δz²)` |
| `triangle-area` | `geometry::triangle_area` | Heron: `16 A² = (a+b+c)(−a+b+c)(a−b+c)(a+b−c)` |
| `corner-angle` | `geometry::corner_angle` | angle `θ ∈ (0, π)` opposite side `a`: `tan²(θ/2) = (a−b+c)(a+b−c) / ((a+b+c)(−a+b+c))` |
| `corner-cot` | `geometry::corner_cotan` | `cot θ = (b² + c² − a²) / (4A)` for the angle opposite side `a` |
| `cotan-laplacian` | `operators::cotan_laplacian` | `L_ij = −½ Σ cot θ` over the corners opposite edge `ij`; `L_ii = −Σ_{j≠i} L_ij` |
| `lumped-mass` | `operators::lumped_mass`, `lumped_mass_normalized` | `m_i = ⅓ Σ A_f` over the faces at `i`; normalized: `m_i / Σ_f A_f`, which sums to 1 |
| `angle-defect` | `operators::angle_defects` | `Ω_i = 2π − Σ θ` over the corners at interior vertex `i` |
| `boundary-geodesic-curvature` | `operators::boundary_geodesic_curvatures` | `k_i = π − Σ θ` over the corners at boundary vertex `i` |
| `gauss-bonnet` | `operators::gauss_bonnet_residual` | `Σ_I Ω_i + Σ_B k_i = 2πχ` (the function returns left minus right) |

## Agreement with Blender's triangulation

`scripts/blender-tessellation-diff.sh [SEED] [SCALE]` lets a headless Blender triangulate random
polygons, then compares its `loop_triangles` with `triangulate_polygons`, triangle by triangle
and corner by corner. Measured with Blender 5.0.1 on Linux x86-64, seeds 1, 2, 3, 11, 12, 13 at
scale 1.0 (87 000 comparisons each):

| polygons | compared | identical |
|----------|---------:|----------:|
| non-degenerate quads (convex, concave, non-planar, random in 3D) | 99 000 | 100 % |
| non-degenerate n-gons, 5 to 400 corners (convex, star-shaped, combs, spirals, non-planar, self-intersecting, random in 3D) | 309 000 | 100 % |
| degenerate quads and n-gons (collinear, coincident or all-collinear corners; rectilinear grid polygons; polygons quantized by single precision) | 114 000 | 100 % |

"Identical" means the same triangles in the same order with the same corner order. Every n-gon
is counted twice, once for each of Blender's two code paths (with and without cached face
normals). A recorded sample of 523 polygons is committed as a test fixture
(`tests/data/blender_tessellation_5.0.1.txt`), so `cargo test` checks the agreement without
Blender.

Known limits of that statement:

- It was measured against one Blender build. The split depends on the signs of single precision
  expressions; a Blender compiled with fused multiply-add contraction (possible on ARM64) can
  decide exact ties differently from this code, which never fuses.
- Blender answers "is a reflex corner inside this ear" with a k-d tree; this crate scans the
  reflex corners. The result is the same; only the running time for huge concave n-gons
  differs. Measured here: a star-shaped polygon with 10 000 corners takes 28 ms, one with
  50 000 corners 0.6 s; a million quads take 90 ms.
- Positions are converted to `f32` for the decision, as Blender stores them. Positions that are
  not exactly representable in `f32` are rounded first.

## Sparse solver timings

`cargo run --release -p ddg-core --example factor_timings` on the largest dataset mesh
(`brain.obj`: 152 893 vertices, 305 782 faces, closed). AMD Threadripper 3970X, faer 0.24.4,
best of 3. `L` is the cotan Laplacian (1.07 M entries), `M` the area-normalized lumped mass, `D`
a diagonal with 301 ones and zeros elsewhere, `K = [[−M, L], [L, D]]` (305 786 unknowns).
"Threads" is the size of the `ThreadPool` passed in (1 = sequential). Times in milliseconds.

| system | method | symbolic | numeric, 1 / 2 / 4 / 8 threads | solve, 1 rhs | residual |
|--------|--------|---------:|-------------------------------:|-------------:|----------|
| `L + 1e-8 I` | `LLᵀ` | 170 | 126 / 111 / 136 / 138 | 20 | 7e-16 |
| `L`, one vertex pinned | `LLᵀ` | 161 | 126 / 113 / 108 / 137 | 20 | 7e-16 |
| `K` | `LDLᵀ`, AMD order | 683 | **fails**: zero pivot | – | – |
| `K` | `LDLᵀ`, paired order | 128 + 98 | 509 / 414 / 331 / 451 | 59 | 1e-8, 5e-16 after 1 refinement step |
| `K` | `LBLᵀ` | 595 | 710 / 1289 / 2090 / 6999 | 71 | **fails**: non-finite solution |
| `K` | `LU` | 863 | 3657 / 2514 / 2573 / 2466 | 97 | 7e-16 |
| `K + diag(0, 1e-8 I)` | `LDLᵀ`, AMD order | 578 | 659 / 530 / 473 / 623 | 67 | 3e-7, 1e-14 after 2 refinement steps |
| `L M⁻¹ L + D` (Schur complement of `K`) | `LLᵀ` | 403 | 559 / 390 / 342 / 448 | 44 | 3e-16 |

Setup on the same mesh: OBJ parse 210 ms, topology 71 ms, geometry 51 ms, Laplacian and mass
74 ms; assembling the Schur complement 335 ms.

What the numbers say:

- Cholesky of the Laplacian costs about 0.3 s once (analysis + factorization) and 0.13 s per
  refactorization; a solve is 20 ms.
- **Threads do not help.** The supernodes of a surface Laplacian are small, and the backend
  parallelizes inside supernodes only. Two to four threads gain at most 15–35 % on the larger
  factors and nothing on the Laplacian; more threads are slower (16 threads: 153 ms for `LLᵀ`
  of `L`; 32 threads: 211–327 ms). `Parallelism::Sequential` is the sensible default.
- **Never hand the backend a large shared pool.** With faer's parallelism pointed at the global
  64-thread rayon pool, `LLᵀ` of `L` took 1.5 s with "16 threads" and 3.4 s with all of them,
  12 to 27 times slower than sequential. This is why the wrapper takes an explicit
  `ThreadPool` and never touches the global pool.
- The saddle-point matrix `K` as given defeats the generic symmetric paths: with the AMD order,
  `LDLᵀ` meets a zero pivot (a zero diagonal entry of `D` is eliminated first), and `LBLᵀ`
  returns a factor whose solve is not finite (its pivoting cannot leave a supernode). `LBLᵀ`
  also gets slower with every added thread (41 s at 32 threads).
- `K` does factorize by plain `LDLᵀ` when unknown `i` of the first block is eliminated
  immediately before unknown `n + i` of the second, with the pairs in the AMD order of `L`
  (`EliminationOrder::Custom`, `fill_reducing_ordering`). No pivot can vanish in that order. Without
  pivoting the solve is accurate to about 1e-8; one step of iterative refinement
  (`Factorization::solve_refined`) restores full accuracy. Cost: 0.23 s once, 0.5 s per
  refactorization, 0.06 s per solve.
- Forming the Schur complement `L M⁻¹ L + D` and using `LLᵀ` costs about the same per
  factorization, but its condition number is that of `L` squared; the residual in the table is
  for the Schur system itself, not for `K`.
- `LU` is the only method that solves `K` to full accuracy with no care at all, at 2.5–3.7 s per
  factorization.

## Tests

```
cargo test -p ddg-core                        # unit, property, dataset and Blender-fixture tests
cargo test -p ddg-core --test dataset -- --nocapture   # prints the dataset topology table
scripts/blender-tessellation-diff.sh 1 1.0    # differential test against an installed Blender
```

The dataset tests read `data/optimal-cones/optimal-cones-data/data` and are skipped, with a
notice, when that directory is absent.

## License

MIT, with one exception: `src/triangulate/blender.rs` is a port of Blender source code
(`mesh_tessellate.cc`, `polyfill_2d.cc` and a few `blenlib` math helpers, Blender v5.0.1) and is
licensed, like its original, under the **GNU General Public License, version 2 or later**
(`LICENSES/GPL-2.0-or-later.txt`). It is compiled only with the cargo feature
`blender-tessellation`, which is enabled by default; a binary built with that feature is a
combined work that must be distributed under the GPL. `--no-default-features` builds the crate
from MIT-licensed code only, without the `triangulate` module.
