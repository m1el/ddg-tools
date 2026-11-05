# ddg-intrinsic-apps

Applications of intrinsic triangulations, Sec. 5 of

> Nicholas Sharp, Yousuf Soliman, Keenan Crane. *Navigating Intrinsic Triangulations.*
> ACM Transactions on Graphics 38(4), 2019. <https://doi.org/10.1145/3306346.3322979>

built on `ddg-intrinsic` (the signpost data structure, Delaunay flipping, Delaunay refinement,
the optimal Delaunay triangulation, edge traces and point queries). The contract is
[`SPEC.md`](SPEC.md) (decisions `N30` to `N49`, continuing those of `ddg-intrinsic`). The paper's
Sec. 5.1 to 5.5 are:

| | module | what |
|---|---|---|
| 5.1 | `steiner` | a short cut through all vertices: the minimal spanning tree of the intrinsic Delaunay triangulation, shortened by Fermat points inserted as intrinsic vertices; its paths on the input mesh; the check that the cut leaves a flat disk |
| 5.2 | `fem` | stiffness, Galerkin and lumped mass matrices on an intrinsic triangulation; the regular 4-1 subdivision; extreme eigenvalues and condition numbers (Lanczos) |
| 5.3 | `heat` | geodesic distance by the heat method on an intrinsic triangulation, factored once, values at input vertices, inserted vertices and any point of the input mesh |
| — | `exact` | the exact polyhedral distance from vertices (window propagation), the ground truth of 5.3; also on a `ddg_intrinsic::Complex` Δ-complex with edge lengths and face-point sources (N37b) |
| 5.4 | `amr` | adaptive refinement: residual error estimator, bulk marking, refinement of the marked faces (`ddg_intrinsic::refine_with`), for the harmonic Green's function, the short-time heat kernel and a harmonic disk map |
| 5.5 | `vector` | tangent fields: transport from the signposts, finite difference and finite element connection Laplacians, harmonic fields with boundary values and the vector maximum principle, smoothest `n`-direction fields and their singularities |

## Methods and papers

- **Steiner trees** (Sec. 5.1). The heuristic of J. M. Smith, D. T. Lee, J. S. Liebman, *An
  O(n log n) heuristic for Steiner minimal tree problems on the Euclidean metric*, Networks
  11(1), 1981, as the paper summarizes it, with the original's priority (the ratio of a
  triangle's Steiner tree to the MST of its corners; described by M. Zachariasen, P. Winter,
  *Concatenation-based greedy heuristics for the Euclidean Steiner tree problem*, DIKU-TR-97/20,
  1997, Algorithmica 25, 1999). The extrinsic comparison is the MST of the input edges, the cut
  of A. Sheffer, J. C. Hart, *Seamster*, IEEE Visualization 2002 (and J. Erickson, S. Har-Peled,
  *Optimally cutting a surface into a disk*, DCG 31, 2004, for the general cut graph). That the
  Delaunay triangulation contains the MST in the plane: G. T. Toussaint, Pattern Recognition 12,
  1980. The Fermat point: E. N. Gilbert, H. O. Pollak, *Steiner minimal trees*, SIAM J. Appl.
  Math. 16, 1968.
- **Finite elements** (Sec. 5.2). The cotan stiffness matrix (R. H. MacNeal, 1949, as the paper
  cites) and the standard Galerkin mass matrix; the condition number by the Lanczos method with
  full reorthogonalization (C. Lanczos 1950; shift-and-invert at zero through a Cholesky
  factorization of the Laplacian with one vertex per component pinned).
- **Geodesic distance** (Sec. 5.3). The heat method: K. Crane, C. Weischedel, M. Wardetzky,
  *Geodesics in heat*, ACM TOG 32(5), 2013, and *The heat method for distance computation*,
  CACM 60(11), 2017. Its source term and the final shift follow geometry-central's
  `heat_method_distance.cpp` (MIT, read, not translated; notice in `LICENSES/`).
- **Exact polyhedral distance** (ground truth). J. S. B. Mitchell, D. M. Mount,
  C. H. Papadimitriou, *The discrete geodesic problem*, SIAM J. Comput. 16(4), 1987;
  V. Surazhsky, T. Surazhsky, D. Kirsanov, S. J. Gortler, H. Hoppe, *Fast exact and approximate
  geodesics on meshes*, ACM TOG 24(3), 2005 (the window representation); S.-Q. Xin, G.-J. Wang,
  *Improving Chen and Han's algorithm on the discrete geodesic problem*, ACM TOG 28(4), 2009 (the
  vertex filter). Written from the papers for this crate.
- **Adaptive refinement** (Sec. 5.4). The residual estimator and the bulk marking of
  P. Morin, R. H. Nochetto, K. G. Siebert, *Convergence of adaptive finite element methods*,
  SIAM Review 44(4), 2002, and K. Mekchay, R. H. Nochetto, *Convergence of adaptive finite
  element methods for general second order linear elliptic PDEs*, SIAM J. Numer. Anal. 43(5),
  2005, as the paper's supplement (Sec. 2) describes them (with `h_T` where the supplement
  prints the area, SPEC N47).
- **Tangent fields** (Sec. 5.5). The finite difference connection Laplacian and the vector
  maximum principle of N. Sharp, Y. Soliman, K. Crane, *The Vector Heat Method*, ACM TOG 38(3),
  2019 (Sec. 5.3, 5.4); the finite element connection Laplacian and the smoothest direction
  fields of F. Knöppel, K. Crane, U. Pinkall, P. Schröder, *Globally Optimal Direction Fields*,
  ACM TOG 32(4), 2013 (Sec. 6, App. D; the local matrices derived anew, SPEC N43).

Downloaded copies of the cited papers that are openly available are in
`../NavigatingIntrinsicTriangulations/cited/applications-5.1-5.3/` and
`../NavigatingIntrinsicTriangulations/cited/applications-5.4-5.5/` (outside the repository).

## Use

```rust
use ddg_core::generators::icosphere;
use ddg_core::Geometry;
use ddg_intrinsic::{DelaunayOptions, IntrinsicTriangulation, RefineOptions};
use ddg_intrinsic_apps::heat::{HeatOptions, HeatSolver, HeatSource};
use ddg_intrinsic_apps::steiner::{SteinerOptions, cut_to_disk, steiner_tree};
use ddg_intrinsic_apps::vector::{DirectionFieldOptions, singularities, smoothest_field};
use ddg_intrinsic_apps::amr::{AmrOptions, Problem, adaptive_solve};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
let sphere = icosphere(3)?;
let mesh = sphere.topology()?;
let geometry = Geometry::from_positions(&mesh, &sphere.positions)?;
let mut tri = IntrinsicTriangulation::new(&mesh, &geometry)?;

// 5.3: flip, refine, factor once, query as often as needed.
tri.flip_to_delaunay(&DelaunayOptions::default())?;
tri.refine(&RefineOptions::default())?;
let solver = HeatSolver::new(&tri, &HeatOptions::default())?;
let d = solver.distance(&[HeatSource::Vertex(0)])?; // d.distances[v], input vertices first

// 5.1: a cut through all vertices that leaves a flat disk.
let tree = steiner_tree(&IntrinsicTriangulation::new(&mesh, &geometry)?, &SteinerOptions::default())?;
let cut = cut_to_disk(&tree.triangulation, &tree.edges)?;
assert!(cut.is_disk && tree.length <= tree.report.mst_length);

// 5.5: the smoothest cross field on the iDT, its singularities (indices add up to 4χ).
let mut idt = IntrinsicTriangulation::new(&mesh, &geometry)?;
idt.flip_to_delaunay(&DelaunayOptions::default())?;
let options = DirectionFieldOptions { n_sym: 4, ..Default::default() };
let (field, report) = smoothest_field(&idt, &options)?;
let sing = singularities(&idt, &field)?;
assert_eq!(sing.index_sum, 8);

// 5.4: the heat kernel of vertex 0, refined where the estimator says.
let mut amr_tri = IntrinsicTriangulation::new(&mesh, &geometry)?;
let problem = Problem::HeatKernel { source: 0, time: 0.01 };
let (solution, amr_report) = adaptive_solve(&mut amr_tri, &problem, &AmrOptions { max_vertices: 5000, ..Default::default() })?;
# Ok(())
# }
```

Python (`ddg_tools`, crate `ddg-py`): `geodesic_distance(positions, triangles, sources, *,
triangulation="refined", time_factor=1.0, boundary="neumann")`, `steiner_tree(positions,
triangles, *, intrinsic=True, steiner_points=True)`, `smoothest_direction_field(positions,
triangles, *, n=1, triangulation="delaunay", laplacian="fem", min_angle_degrees=25.0)` and
`adaptive_solve(positions, triangles, problem, *, source=0, time=None, max_vertices=20000,
bulk=0.5, min_angle_degrees=25.0)`; see `crates/ddg-py/ddg_tools.pyi`.

Blender: nothing. A natural use would be the 5.1 cut as seams for a one-chart unwrap; its
edges are intrinsic (they cross input faces), so seams along them need the cut inserted into
the mesh first (as `ddg_uv::cutting` does for the cut flow). A direction field could drive
Blender's own tools through an attribute. Not done; see the ROADMAP.

## 5.1 Steiner tree approximation

"Many applications require a short cut passing through all vertices of a polyhedron, so that it
can be isometrically unfolded into the plane. A common approach is to cut through extrinsic
edges which approximate the Steiner tree of the vertex set; for a genus-0 surface the shortest
such cut is simply the minimal spanning tree (MST)." The paper flips to the iDT, enqueues all
edges and per-triangle Steiner trees, and adds the shortest element that closes no cycle.
`steiner_tree` does this with the priority of Smith, Lee and Liebman (SPEC N31): a triangle with
two MST sides and all angles below 120° is ranked by `|Steiner tree| / |MST of its corners|`;
all triangles come before the edges, which come in order of length. The Fermat point of every
chosen triangle is inserted (`insert_vertex`), its three edges join the tree. Each such triangle
replaces its two MST sides by a shorter tree, so the result is never longer than the MST of the
same triangulation. The vertices of every boundary loop count as joined (the boundary is cut
already, N31). `cut_to_disk` cuts along the tree and checks: one piece, Euler characteristic 1,
no angle defect off the cut (every curved vertex lies on the cut, so the disk is flat inside),
boundary turning `2π`, and the layout of all triangles by unfolding across the edges that are
not cut agrees with itself.

Lengths on the meshes of the workspace (`examples/apps_study.rs steiner`; lengths in the units
of the mesh). *Extrinsic*: over the input edges, with Fermat points in input faces for the
Steiner column; *intrinsic*: over the iDT. *Flat disk*: `cut_to_disk` found one disk with no
defect inside; bob and the sheared torus have genus one, where a tree does not cut a disk
(Erickson and Har-Peled's cut graphs add loops; not done).

| mesh | V | χ | extrinsic MST | extrinsic Steiner | intrinsic MST | intrinsic Steiner | intrinsic Steiner / extrinsic MST | Fermat points | flat disk |
|---|--:|--:|--:|--:|--:|--:|--:|--:|:-:|
| blub_control_mesh | 112 | 2 | 26.9103 | 26.5116 | 26.8128 | 26.2834 | 0.9767 | 27 | true |
| blub_triangulated | 7106 | 2 | 183.5639 | 182.5806 | 183.3400 | 182.0696 | 0.9919 | 465 | true |
| bob_tri | 5344 | 0 | 144.4965 | 143.8522 | 144.4965 | 143.8974 | 0.9959 | 245 | false |
| bump | 3329 | 1 | 120.0726 | 117.4534 | 120.0726 | 117.4534 | 0.9782 | 428 | true |
| bunny | 14290 | 2 | 236.8960 | 224.0762 | 236.8960 | 224.0762 | 0.9459 | 3898 | true |
| camelbad | 11392 | 1 | 141.3908 | 139.4269 | 139.1484 | 136.4795 | 0.9653 | 2002 | true |
| camelbadnonuniform | 9092 | 1 | 98.4113 | 97.1043 | 96.6217 | 94.7068 | 0.9624 | 1615 | true |
| camel_noise1 | 11381 | 1 | 152.4244 | 149.8478 | 151.5074 | 147.8911 | 0.9703 | 2926 | true |
| camelshotnoise | 11513 | 1 | 87.5276 | 87.0229 | 87.1864 | 86.2599 | 0.9855 | 1358 | true |
| face | 13746 | 1 | 208.3677 | 206.4702 | 208.3597 | 206.0767 | 0.9890 | 1745 | true |
| icosphere_jitter | 642 | 2 | 83.3732 | 78.6090 | 83.3732 | 78.6090 | 0.9429 | 221 | true |
| sheared_square | 289 | 1 | 14.0625 | 14.0625 | 14.0625 | 14.0625 | 1.0000 | 0 | true |
| sheared_torus | 960 | 0 | 106.8226 | 106.8226 | 102.0730 | 101.8591 | 0.9535 | 15 | false |
| spot_control_mesh | 188 | 2 | 29.1235 | 28.5730 | 29.0954 | 28.4096 | 0.9755 | 48 | true |
| spot_triangulated | 2930 | 2 | 103.4547 | 102.0594 | 103.4547 | 102.0081 | 0.9860 | 363 | true |


The intrinsic Steiner tree is 1 to 6 % shorter than the extrinsic MST (spot 1.4 %, camelbad
3.5 %, the jittered icosphere 5.7 %, the bunny 5.4 %). The Fermat points do most of it; the
intrinsic MST is shorter than the extrinsic one only where the iDT has edges that the input
lacks (camelbad 1.6 %, the sheared torus 4.4 %), and equal on meshes that are Delaunay already
(bump, bunny, the icosphere). The paper reports no table; its Fig. 15 shows a cut that is
shorter than the MST, and it states that the ratio approaches 20 % on a constructed example
(moving two vertices of its Fig. 15 apart). In the plane Smith, Lee and Liebman's heuristic
saves about 3 % over the MST on random points (Zachariasen and Winter, Sec. 1), which is the
range measured here. On the sheared square the boundary is joined from the start and the MST
becomes a set of combs of lattice edges to the boundary; no triangle has two MST sides and no
Fermat point is inserted.

## 5.2 Finite elements

"The stiffness matrix is built exactly as with an ordinary triangle mesh, via the cotan
formula; we likewise use the standard Galerkin mass matrix. For an equal number of vertices,
the intrinsic retriangulation (via either iDR or iODT) yields better condition numbers
(Figure 17) than an extrinsic version of ODT [...]; values for regular 4-1 refinement are off
the chart." `stiffness_matrix`, `galerkin_mass_matrix`, `lumped_mass` assemble from the
intrinsic edge lengths; `regular_subdivision` is the extrinsic 4-1 refinement, from lengths
alone (halves of the sides, midlines half the opposite side: the same surface and angles);
`condition_number` gives `λ_max / λ_min` of `L x = λ M x` over the complement of the constants,
for `M` the Galerkin mass, the lumped mass or the identity (SPEC N34). The paper does not say
which matrix its Fig. 17 conditions; all three are reported.

Coarse inputs with bad triangles (`examples/apps_study.rs fem`): the input, its iDT, and at the
vertex counts of one and two levels of 4-1 subdivision: the subdivision, iDR (25°, with the
circumradius bound found by bisection to match the count) and iODT (target length found by
bisection). The extrinsic ODT of the figure is not implemented.

| mesh | triangulation | V | min angle | max angle | κ(L, M_Galerkin) | κ(L, M_lumped) | κ(L) |
|---|---|--:|--:|--:|--:|--:|--:|
| blub_control_mesh | input | 112 | 2.09° | 175.04° | 1.873e4 | 7.527e3 | 1.174e3 |
| blub_control_mesh | idt | 112 | 3.82° | 145.83° | 1.026e4 | 4.525e3 | 7.212e2 |
| blub_control_mesh | 4-1 x1 | 442 | 2.09° | 175.04° | 1.669e5 | 6.907e4 | 8.583e3 |
| blub_control_mesh | iDR @4-1 x1 | 440 | 25.21° | 129.58° | 1.892e4 | 7.880e3 | 6.184e2 |
| blub_control_mesh | iODT @4-1 x1 | 443 | 12.13° | 130.64° | 8.437e3 | 3.582e3 | 1.114e3 |
| blub_control_mesh | 4-1 x2 | 1762 | 2.09° | 175.04° | 1.367e6 | 4.900e5 | 4.774e4 |
| blub_control_mesh | iDR @4-1 x2 | 1735 | 25.21° | 129.58° | 1.811e4 | 6.145e3 | 2.342e3 |
| blub_control_mesh | iODT @4-1 x2 | 1772 | 22.41° | 110.37° | 1.295e4 | 4.824e3 | 2.464e3 |
| bob_controlmesh | input | 334 | 7.94° | 140.72° | 5.986e3 | 2.207e3 | 3.376e2 |
| bob_controlmesh | idt | 334 | 8.41° | 110.45° | 1.013e4 | 4.145e3 | 3.246e2 |
| bob_controlmesh | 4-1 x1 | 1336 | 7.94° | 140.72° | 4.476e4 | 1.614e4 | 1.727e3 |
| bob_controlmesh | iDR @4-1 x1 | 1331 | 25.04° | 113.14° | 1.543e4 | 6.312e3 | 7.196e2 |
| bob_controlmesh | iODT @4-1 x1 | 1348 | 20.52° | 112.02° | 1.221e4 | 4.697e3 | 8.547e2 |
| bob_controlmesh | 4-1 x2 | 5344 | 7.94° | 140.72° | 2.423e5 | 8.004e4 | 8.112e3 |
| bob_controlmesh | iDR @4-1 x2 | 5300 | 30.98° | 116.10° | 1.642e4 | 6.308e3 | 3.190e3 |
| bob_controlmesh | iODT @4-1 x2 | 5412 | 23.89° | 109.37° | 2.235e4 | 7.577e3 | 3.313e3 |
| bump | input | 3329 | 25.37° | 95.94° | 1.783e4 | 5.913e3 | 1.993e3 |
| bump | idt | 3329 | 25.37° | 95.94° | 1.783e4 | 5.913e3 | 1.993e3 |
| bump | 4-1 x1 | 13185 | 25.37° | 95.94° | 8.372e4 | 2.661e4 | 8.448e3 |
| bump | iDR @4-1 x1 | 13118 | 25.66° | 118.97° | 5.632e4 | 2.238e4 | 8.283e3 |
| bump | iODT @4-1 x1 | 13157 | 26.89° | 112.71° | 4.945e4 | 1.740e4 | 7.257e3 |
| bump | 4-1 x2 | 52481 | 25.37° | 95.94° | 3.721e5 | 1.119e5 | 3.556e4 |
| bump | iDR @4-1 x2 | 52218 | 25.56° | 119.30° | 1.775e5 | 6.836e4 | 3.469e4 |
| bump | iODT @4-1 x2 | 51521 | 23.10° | 113.77° | 2.186e5 | 6.987e4 | 3.150e4 |
| icosphere_jitter | input | 642 | 23.45° | 117.30° | 6.395e2 | 1.994e2 | 1.808e2 |
| icosphere_jitter | idt | 642 | 23.45° | 117.30° | 6.395e2 | 1.994e2 | 1.808e2 |
| icosphere_jitter | 4-1 x1 | 2562 | 23.45° | 117.30° | 3.623e3 | 1.384e3 | 9.271e2 |
| icosphere_jitter | iDR @4-1 x1 | 2566 | 30.50° | 117.75° | 3.621e3 | 1.177e3 | 8.320e2 |
| icosphere_jitter | iODT @4-1 x1 | 2578 | 31.29° | 114.76° | 3.541e3 | 1.263e3 | 6.773e2 |
| icosphere_jitter | 4-1 x2 | 10242 | 23.45° | 117.30° | 2.425e4 | 8.808e3 | 4.700e3 |
| icosphere_jitter | iDR @4-1 x2 | 10368 | 30.23° | 118.49° | 1.436e4 | 4.969e3 | 3.314e3 |
| icosphere_jitter | iODT @4-1 x2 | 10188 | 27.27° | 104.57° | 1.357e4 | 4.732e3 | 2.774e3 |
| sheared_square | input | 289 | 2.73° | 165.96° | 1.506e5 | 4.580e4 | 5.277e4 |
| sheared_square | idt | 289 | 4.40° | 161.57° | 9.521e4 | 2.709e4 | 1.110e4 |
| sheared_square | 4-1 x1 | 1089 | 2.73° | 165.96° | 6.298e5 | 1.903e5 | 1.987e5 |
| sheared_square | iDR @4-1 x1 | 1165 | 14.04° | 102.68° | 6.018e5 | 1.710e5 | 1.643e4 |
| sheared_square | iODT @4-1 x1 | 865 | 14.04° | 151.93° | 2.418e5 | 7.234e4 | 1.349e4 |
| sheared_square | 4-1 x2 | 4225 | 2.73° | 165.96° | 2.572e6 | 7.758e5 | 7.721e5 |
| sheared_square | iDR @4-1 x2 | 4205 | 14.04° | 112.49° | 2.405e6 | 6.840e5 | 5.903e4 |
| sheared_square | iODT @4-1 x2 | 4341 | 14.04° | 118.15° | 1.751e6 | 5.321e5 | 6.340e4 |
| sheared_torus | input | 960 | 0.56° | 176.33° | 3.261e4 | 1.211e4 | 6.779e3 |
| sheared_torus | idt | 960 | 16.54° | 91.40° | 4.100e3 | 1.335e3 | 5.743e2 |
| sheared_torus | 4-1 x1 | 3840 | 0.56° | 176.33° | 1.824e5 | 6.131e4 | 3.601e4 |
| sheared_torus | iDR @4-1 x1 | 3869 | 26.65° | 111.96° | 5.925e3 | 1.901e3 | 1.540e3 |
| sheared_torus | iODT @4-1 x1 | 3900 | 16.09° | 113.25° | 1.860e4 | 5.931e3 | 2.504e3 |
| sheared_torus | 4-1 x2 | 15360 | 0.56° | 176.33° | 4.073e6 | 1.626e6 | 2.362e5 |
| sheared_torus | iDR @4-1 x2 | 15120 | 31.10° | 117.80° | 2.525e4 | 7.985e3 | 5.182e3 |
| sheared_torus | iODT @4-1 x2 | 15180 | 31.67° | 93.53° | 2.162e4 | 7.071e3 | 5.240e3 |
| spot_control_mesh | input | 188 | 0.17° | 173.32° | 7.193e4 | 3.106e4 | 9.039e3 |
| spot_control_mesh | idt | 188 | 1.49° | 140.61° | 1.468e4 | 5.132e3 | 1.638e3 |
| spot_control_mesh | 4-1 x1 | 746 | 0.17° | 173.32° | 1.606e7 | 7.681e6 | 6.956e4 |
| spot_control_mesh | iDR @4-1 x1 | 739 | 25.01° | 121.74° | 3.585e5 | 1.422e5 | 9.416e2 |
| spot_control_mesh | iODT @4-1 x1 | 735 | 1.49° | 123.98° | 2.437e4 | 1.019e4 | 5.941e3 |
| spot_control_mesh | 4-1 x2 | 2978 | 0.17° | 173.32° | 6.682e8 | 2.747e8 | 4.819e5 |
| spot_control_mesh | iDR @4-1 x2 | 2936 | 25.25° | 123.19° | 3.130e5 | 1.149e5 | 2.971e3 |
| spot_control_mesh | iODT @4-1 x2 | 2994 | 2.86° | 137.54° | 4.490e4 | 1.920e4 | 1.326e4 |


The 4-1 subdivision keeps the angles of the input, and its condition number grows with the
vertex count from the input's on: on spot's control mesh (corners of 0.17°) `κ(L, M_G)` is
`6.7e8` at 2978 vertices, against `3.1e5` (iDR) and `4.5e4` (iODT) at about the same count; on
the sheared torus (176°) `4.1e6` against `2.5e4` and `2.2e4`; on blub's control mesh `1.4e6`
against `1.8e4` and `1.3e4`. That is the figure's finding ("off the chart"). Where the input
is well shaped (bump, the jittered icosphere) the three are within a factor of two. Between iDR
and iODT neither is better everywhere: iODT is lower in 11 of the 14 matched comparisons of the
Galerkin column, iDR in 3 (bob and bump at two levels, the sheared torus at one level); on the
sheared square the boundary keeps both above `1e5` (its input corners of 14° stay). The match
is loose in one row: on the sheared square at one level the bisection of iODT's target length
missed the 1089 vertices of the subdivision and ended at 865 (iDR has 1165).

## 5.3 Geodesic distance

"PDE-based methods for computing geodesic distance such as the heat method are attractive
because they are easy to implement and can take advantage of fast linear solvers [...] by just
using the iDT the error is dramatically reduced; further refinement with iDR reduces the error
to less than 1%." (Fig. 18: mean error 59.6 % on the input, 20.4 % on the iDT, 0.7 % after iDR.)

`HeatSolver::new(tri, &HeatOptions)` assembles and factors `M + t L` (lumped mass, `t = m h²`,
`h` the mean edge length of the triangulation, `m = 1`) and the Laplacian with one vertex per
component pinned; `distance(&[HeatSource])` solves the heat step, normalizes `−∇u` per triangle,
integrates its divergence and solves the Poisson step, then shifts the result to vanish at the
sources (SPEC N35). Sources are vertices of the triangulation or points of the input mesh
(found by a point query). `value_at(values, point)` evaluates the result at any point of the
input mesh; inserted vertices have their place on the input mesh in `vertex_location` (N36).

Errors against the exact polyhedral distance (`exact::exact_distance`, N37) at the input
vertices, three sources per mesh (vertices `0`, `n/3`, `2n/3`; `examples/apps_study.rs heat`):
the mean of `|d − d*| / d*`, the mean of `|d − d*|` over the largest distance, and the largest
`|d − d*|` over the largest distance. *iDR*: refinement to 25°; *iDR r/2*: also every
circumradius at most half the mean edge of the iDT (N38). Times: factorization (Laplacian, mass,
both Cholesky factorizations) and one query, single thread, one process at a time, on this
machine (AMD, 64 cores, of which one was used).

| mesh | triangulation | V | factor (s) | query (s) | mean relative error | mean error / max d | max error / max d |
|---|---|--:|--:|--:|--:|--:|--:|
| blub_triangulated | input | 7106 | 0.0299 | 0.00271 | 1.36 % | 0.70 % | 4.58 % |
| blub_triangulated | iDT | 7106 | 0.0274 | 0.00264 | 1.28 % | 0.63 % | 3.22 % |
| blub_triangulated | iDR | 8877 | 0.0334 | 0.00335 | 1.73 % | 0.89 % | 3.88 % |
| blub_triangulated | iDR r/2 | 17510 | 0.1008 | 0.00891 | 0.86 % | 0.41 % | 1.93 % |
| bob_tri | input | 5344 | 0.0519 | 0.00228 | 2.28 % | 1.17 % | 4.11 % |
| bob_tri | iDT | 5344 | 0.0214 | 0.00208 | 2.26 % | 1.16 % | 4.07 % |
| bob_tri | iDR | 5437 | 0.0216 | 0.00207 | 2.24 % | 1.15 % | 3.98 % |
| bob_tri | iDR r/2 | 10796 | 0.0381 | 0.00547 | 1.41 % | 0.69 % | 2.98 % |
| bump | input | 3329 | 0.0124 | 0.00120 | 1.12 % | 0.47 % | 2.15 % |
| bump | iDT | 3329 | 0.0114 | 0.00119 | 1.12 % | 0.47 % | 2.15 % |
| bump | iDR | 3329 | 0.0116 | 0.00121 | 1.12 % | 0.47 % | 2.15 % |
| bump | iDR r/2 | 8210 | 0.0422 | 0.00333 | 0.71 % | 0.28 % | 2.00 % |
| camelbad | input | 11392 | 0.0499 | 0.00432 | 14.61 % | 8.45 % | 32.19 % |
| camelbad | iDT | 11392 | 0.0445 | 0.00439 | 1.59 % | 0.80 % | 4.13 % |
| camelbad | iDR | 15561 | 0.0611 | 0.00603 | 1.56 % | 0.78 % | 5.65 % |
| camelbad | iDR r/2 | 30595 | 0.1582 | 0.01618 | 0.80 % | 0.39 % | 2.54 % |
| camelbadnonuniform | input | 9092 | 0.0345 | 0.00339 | 22.44 % | 11.20 % | 74.26 % |
| camelbadnonuniform | iDT | 9092 | 0.0305 | 0.00329 | 4.18 % | 1.84 % | 30.62 % |
| camelbadnonuniform | iDR | 12848 | 0.0464 | 0.00494 | 2.93 % | 1.25 % | 27.24 % |
| camelbadnonuniform | iDR r/2 | 28696 | 0.1460 | 0.01515 | 0.85 % | 0.31 % | 2.09 % |
| camel_noise1 | input | 11381 | 0.0708 | 0.00548 | 2.18 % | 1.02 % | 6.82 % |
| camel_noise1 | iDT | 11381 | 0.0516 | 0.00432 | 2.07 % | 1.04 % | 4.58 % |
| camel_noise1 | iDR | 18079 | 0.0932 | 0.00714 | 2.23 % | 1.10 % | 5.65 % |
| camel_noise1 | iDR r/2 | 30492 | 0.1428 | 0.01552 | 1.58 % | 0.81 % | 3.58 % |
| face | input | 13746 | 0.0860 | 0.00620 | 1.00 % | 0.45 % | 2.27 % |
| face | iDT | 13746 | 0.0488 | 0.00657 | 0.77 % | 0.36 % | 2.04 % |
| face | iDR | 13850 | 0.0493 | 0.00659 | 0.78 % | 0.36 % | 2.04 % |
| face | iDR r/2 | 28436 | 0.1190 | 0.01518 | 0.65 % | 0.30 % | 1.69 % |
| icosphere_jitter | input | 642 | 0.0020 | 0.00023 | 1.54 % | 0.56 % | 1.31 % |
| icosphere_jitter | iDT | 642 | 0.0022 | 0.00023 | 1.54 % | 0.56 % | 1.30 % |
| icosphere_jitter | iDR | 643 | 0.0019 | 0.00022 | 1.54 % | 0.56 % | 1.30 % |
| icosphere_jitter | iDR r/2 | 1572 | 0.0061 | 0.00056 | 1.20 % | 0.48 % | 1.29 % |
| sheared_square | input | 289 | 0.0006 | 0.00008 | 60.46 % | 30.29 % | 109.05 % |
| sheared_square | iDT | 289 | 0.0004 | 0.00008 | 52.74 % | 27.96 % | 85.37 % |
| sheared_square | iDR | 543 | 0.0008 | 0.00013 | 1.42 % | 0.42 % | 1.16 % |
| sheared_square | iDR r/2 | 705 | 0.0011 | 0.00016 | 1.10 % | 0.28 % | 0.89 % |
| sheared_torus | input | 960 | 0.0032 | 0.00033 | 27.55 % | 13.86 % | 52.81 % |
| sheared_torus | iDT | 960 | 0.0025 | 0.00032 | 4.04 % | 1.97 % | 9.29 % |
| sheared_torus | iDR | 1221 | 0.0035 | 0.00041 | 4.04 % | 1.88 % | 9.29 % |
| sheared_torus | iDR r/2 | 1916 | 0.0075 | 0.00069 | 2.59 % | 1.26 % | 5.72 % |
| spot_control_mesh | input | 188 | 0.0005 | 0.00006 | 5.84 % | 3.05 % | 10.44 % |
| spot_control_mesh | iDT | 188 | 0.0004 | 0.00006 | 5.69 % | 2.90 % | 10.12 % |
| spot_control_mesh | iDR | 533 | 0.0014 | 0.00017 | 4.41 % | 2.15 % | 8.29 % |
| spot_control_mesh | iDR r/2 | 715 | 0.0018 | 0.00023 | 3.58 % | 1.78 % | 8.72 % |
| spot_triangulated | input | 2930 | 0.0105 | 0.00103 | 1.80 % | 0.95 % | 4.58 % |
| spot_triangulated | iDT | 2930 | 0.0094 | 0.00102 | 1.78 % | 0.94 % | 4.64 % |
| spot_triangulated | iDR | 3185 | 0.0100 | 0.00110 | 1.80 % | 0.95 % | 4.55 % |
| spot_triangulated | iDR r/2 | 7073 | 0.0302 | 0.00277 | 1.11 % | 0.61 % | 4.54 % |

Where the input has bad triangles the paper's finding reproduces: camelbad 14.6 % on the input,
1.6 % on the iDT, 0.8 % after iDR with the size bound; camelbadnonuniform 22.4 %, 4.2 %, 0.85 %;
the sheared torus 27.6 %, 4.0 %, 2.6 %. On well-shaped meshes (spot, bob, blub, face, the
icosphere) the input is already Delaunay or nearly so and all four are within the heat method's
own error (1 to 2.3 %), which shrinks only with the mesh width: refinement to an angle bound
alone does not change it (it adds few vertices; cf. `ddg-intrinsic` SPEC N16), the size bound
does. The paper's "less than 1 %" after iDR is reached here with the size bound (camelbad,
camelbadnonuniform, face, blub), not with the angle bound alone.

*Surfaces with boundary.* On the sheared square (a flat parallelogram of obtuse triangles)
the iDT barely helps (60.5 % → 52.7 %): boundary edges cannot be flipped, the iDT keeps a row
of slivers (angles up to 161°) along the two slanted sides, and in a sliver the gradient of the
heat is dominated by the error across it; the divergence carries that error into the whole
domain. Refinement splits the encroached boundary edges (`ddg-intrinsic` SPEC N1) and brings the
error to 1.4 %. The time factor `m` matters there too, in a check made at the far corner of
the square only (not a table): with `m = 10` the iDT's distance there was right (a smoother
field), with `m = 0.1` far worse.

The bench `benches/heat.rs` (criterion) measures the factorization and one query on an
icosphere, spot and camelbad, on the input, the iDT and iDR. Run alone on a quiet machine (AMD Ryzen Threadripper 3970X (32 cores, 64 threads), 2026-10-19, one process at a time pinned to one core, load average 1.0 to 1.5;
`runs_tmp/bench-quiet/heat.log`), factorization / one distance query (median):

| mesh | input | iDT | iDR (25°) |
|---|---|---|---|
| icosphere4 | 8.4 ms / 0.87 ms | 9.5 ms / 0.88 ms | 8.7 ms / 0.87 ms |
| spot | 6.9 ms / 0.93 ms | 7.1 ms / 0.95 ms | 7.7 ms / 1.03 ms |
| camelbad | 46.6 ms / 3.89 ms | 46.7 ms / 3.90 ms | 62.0 ms / 5.45 ms |

(The table above was also measured one process at a time.)

## 5.4 Adaptive mesh refinement

"The basic idea is to use a local a posteriori error estimate as an additional refinement
criterion in iDR [...] adaptive refinement is 2–10 times faster than uniform refinement (taking
all computation into account)" (Fig. 19: for the same accuracy uniform iDR needs 18× the
vertices for the harmonic Green's function and 54× for the short-time heat kernel). Fig. 20: the
harmonic map of a disk flips triangles on a bad input, not on the iDT, and AMR "intelligently
reduces distortion".

`adaptive_solve(tri, &problem, &AmrOptions)` refines to the angle bound, then repeats: solve
(`solve`: cotan stiffness, Galerkin mass), estimate (`estimate`: the residual indicator `η_T`
and the oscillation of Mekchay and Nochetto, with `h_T` the diameter; SPEC N47), mark (`mark`:
the supplement's bulk criterion with `θ = 1/2`, or Dörfler's on squares; N48), refine the marked
faces with `IntrinsicTriangulation::refine_with` (N48). The problems (N46): the Green's function
`Δu = δ_s − 1/|M|`, one backward Euler step `(1 + tΔ)u = δ_s` with `t` the squared mean input
edge length, and a harmonic map onto the unit disk.

Measured by `examples/amr_study.rs` (N49): relative `L²` error against a reference of about
700 000 vertices (adaptive, then uniformly refined; a half-size reference differs from it by
`7e-7` to `5e-4` relative), source at vertex 0 of the control meshes. For one step of the
adaptive run near the paper's accuracy (`7e-3`), the uniform iDR with the same error
(interpolated log–log between uniform runs):

| mesh | problem | marking | error | AMR V | AMR s | uniform V | uniform s | V ratio | time ratio |
|---|---|---|--:|--:|--:|--:|--:|--:|--:|
| spot_control_mesh | Green | supplement | 6.87e-3 | 3533 | 0.066 | 11648 | 0.072 | 3.3 | 1.1 |
| bob_controlmesh | Green | supplement | 6.66e-3 | 2440 | 0.044 | 45458 | 0.299 | 18.6 | 6.7 |
| blub_control_mesh | Green | supplement | 6.48e-3 | 5088 | 0.082 | 30739 | 0.186 | 6.0 | 2.3 |
| spot_control_mesh | Green | squares | 7.18e-3 | 5814 | 0.382 | 10503 | 0.066 | 1.8 | 0.2 |
| spot_control_mesh | heat | supplement | 7.08e-3 | 1740 | 0.037 | 152713 | 1.579 | 87.8 | 43.1 |
| bob_controlmesh | heat | supplement | 1.09e-2 | 1569 | 0.028 | 94375 | 0.847 | 60.2 | 29.8 |
| blub_control_mesh | heat | supplement | 5.87e-3 | 2624 | 0.044 | 134733 | 1.309 | 51.3 | 29.9 |
| spot_control_mesh | heat | squares | 7.16e-3 | 1758 | 0.086 | 147999 | 1.556 | 84.2 | 18.0 |

(bob's heat kernel: the last adaptive step within the range of the uniform runs. *AMR s*: all
steps of the adaptive run including the initial refinement; *uniform s*: one uniform run,
refinement and solve. The times were re-measured alone on a quiet machine (AMD Ryzen Threadripper 3970X (32 cores, 64 threads), 2026-10-19, one process at a time pinned to one core, load average 1.0 to 1.5; single thread, `runs_tmp/bench-quiet/amr-*.log`); the vertex
counts and errors reproduced exactly.)

The heat kernel reproduces the paper (51× to 88× the vertices, against its 54×); the Green's
function needs 3× to 19× (the paper's 18× is the upper end), and its time ratio spans 1.1× to
6.7× (paper: 2–10×). The supplement's criterion (sums of indicators) beats Dörfler's on squares
clearly for the Green's function: the squares mark few faces per step and need many more steps
(the spot rows: 0.38 s against 0.07 s). The study prints the full table per step.

Disk maps (Fig. 20), `jittered_disk` (a jittered sheared grid, 121 vertices) and `bumpy_disk`
(a sheared height field, 169 vertices): flipped faces of the map and the conformal distortion
`σ₁/σ₂` (area-weighted mean and 99 % quantile), AMR at its last step:

| mesh | triangulation | V | L² error | flipped | mean σ₁/σ₂ | q99 σ₁/σ₂ |
|---|---|--:|--:|--:|--:|--:|
| jittered_disk | input | 121 | 4.6e-2 | 1 | 8.17 | 34.6 |
| jittered_disk | iDT | 121 | 4.4e-3 | 0 | 8.39 | 37.8 |
| jittered_disk | AMR | 20296 | 1.1e-5 | 0 | 8.38 | 35.3 |
| bumpy_disk | input | 169 | 1.9e-2 | 0 | 2.55 | 12.7 |
| bumpy_disk | iDT | 169 | 1.1e-2 | 0 | 2.59 | 18.1 |
| bumpy_disk | AMR | 17435 | 8.2e-5 | 0 | 2.44 | 7.9 |

The input of the jittered disk flips a face, the iDT and every AMR step do not (injective, as
the paper says: non-negative weights and a convex boundary polygon). AMR lowers the 99 %
quantile of the distortion on the bumpy disk (12.7 → 7.9) but not on the jittered one, where
the distortion is that of mapping a long parallelogram onto a disk. For the same error AMR
needs 1.1× to 1.6× fewer vertices than uniform refinement and is not faster here (time ratio
0.5–0.8). The largest distortion grows with refinement (to about 200–490) at the corners of the
parallelogram, where the map is singular.

## 5.5 Tangent vector fields

"Since our signpost data structure already maintains tangent coordinate systems across all
operations, it is naturally suited for problems involving tangent vector fields." Fig. 21: the
smoothest vector field with the finite element connection Laplacian of Knöppel et al. has
Dirichlet energy 1200.2 on the input, 123.5 on the iDT, 47.0 on the iDR. Fig. 22: with
prescribed boundary vectors, the harmonic field with the finite difference connection
Laplacian satisfies the vector maximum principle on the iDT (every vector in the convex cone of
its transported neighbours) and may not on the input.

`connection_matrices` builds both Laplacians from the signposts, or from the same
input-aligned frames computed on integer coordinates (`vector::TangentFrames`, below; transport
`ρ_ij = φ_ji + π − φ_ij`, SPEC N40 to N43), `smoothest_field` the smallest eigenvector by block
inverse iteration (N44), `singularities` the indices, `harmonic_field` and
`check_maximum_principle` the rest (N45).

Smoothest fields (`examples/vector_study.rs smoothest`): energy `λ · area` (finite elements;
the field normalized to mean square length one) on the input, the iDT and the iDR (25°):

| mesh | V input / iDR | n | input | iDT | iDR | singular faces input / iDT / iDR |
|---|--:|--:|--:|--:|--:|--:|
| spot_control_mesh | 188 / 533 | 1 | 44.27 | 42.25 | 42.02 | 6 / 10 / 14 |
| spot_control_mesh | | 4 | 100.23 | 96.00 | 96.34 | 54 / 46 / 54 |
| bob_controlmesh | 334 / 374 | 1 | 17.66 | 17.48 | 17.58 | 2 / 2 / 2 |
| spot_triangulated | 2930 / 3185 | 1 | 30.45 | 30.42 | 30.43 | 8 / 8 / 8 |
| bob_tri | 5344 / 5437 | 1 | 16.76 | 16.76 | 16.77 | 2 / 2 / 2 |
| blub_triangulated | 7106 / 8877 | 1 | 30.60 | 30.59 | 30.63 | 8 / 8 / 8 |
| bunny | 14290 / 14290 | 1 | 25.49 | 25.49 | 25.49 | 4 / 4 / 4 |
| cow1 | 15049 / 15627 | 1 | 14.83 | 14.83 | 14.83 | 9 / 9 / 9 |
| airplane | 2349 / 3558 | 1 | 24.66 | 24.53 | 24.29 | 9 / 9 / 9 |
| camelbad | 11392 / 15561 | 1 | 20.99 | 20.77 | 20.78 | 1 / 1 / 1 |
| camelbadnonuniform | 9092 / 12848 | 1 | 21.39 | 20.92 | 20.61 | 1 / 1 / 3 |
| camelshotnoise | 11513 / 17173 | 1 | 60.85 | 58.59 | 34.55 | 111 / 100 / 39 |
| camelshotnoise | | 4 | 86.26 | 82.93 | 50.66 | 143 / 188 / 141 |
| camel_noise2 | 11381 / 19749 | 1 | 154.28 | 146.29 | 159.70 | 13 / 17 / 15 |

On a closed surface the indices add up to `nχ` on every triangulation (spot 2, 4, 8 for
`n = 1, 2, 4`; bob 0); on a surface with boundary the sum is not a topological invariant and
changes with the triangulation. The block solver converged on every run (16 to 64 iterations).
A round sphere gives `λ ≈ n` as in the smooth case (`tests/vector.rs`). Times per field, alone on a quiet machine
(AMD Ryzen Threadripper 3970X (32 cores, 64 threads), 2026-10-19, one process at a time pinned to one core, load average 1.0 to 1.5; `runs_tmp/bench-quiet/vector.log`): 0.013 s (spot's control mesh) to 2.7 s (the bunny's cross field); per mesh, the
range over the 1-, 2- and 4-direction fields on the three triangulations: spot 0.23 to 0.29 s,
bob_tri 0.36 to 0.90 s, blub 0.58 to 0.91 s, camelbad 0.89 to 1.5 s, camel_noise2 0.89 to 2.5 s.

Fig. 21's tenfold drop does not appear on these meshes: where the input is well shaped (spot,
bob, blub, bunny, cow) the three energies agree within 1 %, and on the bad camels the iDT
lowers the energy by 1–5 %. The one large drop is on the shot-noise camel at iDR (`n = 1`:
60.9 → 58.6 → 34.6, 43 % below the input), whose spikes the refinement resolves; on
camel_noise2 iDR raises the energy above the input (146 → 160) while adding 70 % vertices. The
mesh of Fig. 21 is not in the data set, so its numbers were not reproduced; the paper itself
notes that the iDT need not minimize the vector Dirichlet energy.

Harmonic fields with the inward normal on the boundary (`examples/vector_study.rs harmonic`,
finite differences as in the paper; `sheared:n:shear:jitter:seed` is an `n × n` grid of the
unit square, jittered and sheared): negative cotan weights of the triangulation, interior
vertices checked, vertices outside the cone of their neighbours, and the smallest margin
(radians; N45):

| mesh | triangulation | V | negative weights | checked | outside the cone | min margin |
|---|---|--:|--:|--:|--:|--:|
| sheared 12, shear 1.5 | input | 169 | 133 | 121 | 0 | 0.719 |
| | iDT | 169 | 0 | 121 | 0 | 0.417 |
| sheared 12, shear 4 | input | 169 | 131 | 121 | 10 | −0.120 |
| | iDT | 169 | 0 | 121 | 0 | 0.000 |
| | iDR | 454 | 1 | 257 | 0 | 0.000 |
| sheared 16, shear 6 | input | 289 | 248 | 225 | 11 | −0.040 |
| | iDT | 289 | 0 | 225 | 0 | 0.000 |
| | iDR | 966 | 5 | 579 | 0 | 0.000 |
| sheared 20, shear 8 | input | 441 | 361 | 361 | 12 | −0.016 |
| | iDT | 441 | 0 | 361 | 0 | 0.000 |
| | iDR | 1613 | 8 | 964 | 0 | 0.000 |
| bump | input = iDT = iDR | 3329 | 0 | 3201 | 0 | 2.559 |
| face | input | 13746 | 6539 | 13406 | 0 | 1.996 |
| | iDT | 13746 | 0 | 13406 | 0 | 2.344 |
| | iDR | 13850 | 0 | 13466 | 0 | 2.298 |

The maximum principle fails on the strongly sheared inputs (10–12 vertices; 8–17 with the
finite element Laplacian) and holds on every iDT and iDR, as Fig. 22 shows. A margin of `0.000`
is a vector parallel to a neighbour's; the few "negative" weights of the iDR are `−1e-15` and
smaller, rounding. It holds on the face's input too, despite 6539 negative weights: the
condition is sufficient, not necessary.

## Representations (Integer Coordinates, WP11)

`fem`, `heat`, `steiner` and `amr` are written over `ddg_intrinsic::Triangulation` and run on
signposts (`IntrinsicTriangulation`) and on integer coordinates
(`ddg_integer_coords::IntegerTriangulation`, Gillespie, Sharp, Crane 2021). The matrices come
from `ddg_intrinsic::generic_ops::{cotan_laplacian, lumped_mass}` (bitwise the signposts'
methods); point sources of the heat method and the cross-representation error of `amr` need
`representation::InputPoints` (`T^0 -> T^1`). `representation::Representation::DEFAULT` is
**`Integer`**, since the Thingi10k evaluation (WP10) passed on all 7745 meshes (Integer
Coordinates SPEC §13 Q3); the Python and Blender front ends follow it. Signposts stay an option
(`Representation::Signposts`). This crate depends on `ddg-integer-coords`; the reverse dependency must
not exist.

`vector` (tangent fields and connection Laplacians) is generic over `vector::TangentFrames`,
the direction of every halfedge at its tail in the input-aligned frame of the tail: signposts
borrow their stored angles (every bit as before, pinned by `tests/vector_pins.rs`); integer
coordinates compute them from the roundabouts and the straight segments of the edges over the
input mesh (`IntegerTriangulation::frame_directions`, Integer Coordinates SPEC §4.12, I26;
`AngleFrames` are not input-aligned). The functions that cannot fail on signposts keep their
signatures (`transport_angle`, `halfedge_directions`, `check_maximum_principle`) next to
generic ones that can (`transport_angles`, `halfedge_directions_of`,
`check_maximum_principle_of`). `tests/vector_representation.rs` (release; jittered icospheres,
bumpy and sheared grids, sheared torus, polar disk, face, camelbad, spot; iDT and Ruppert
refinement on both, 22 pairs with the same triangulation): directions and transport within
`3.5e-11` rad (spot refined; `1e-13` or less on the others), stiffness matrices within `2.8e-11`
relative, smoothest fields (`n` = 1, 2, 4, both Laplacians, 68 runs) within `5.7e-10` after
aligning the phase with the same eigenvalues (relative `1e-14`) and **the same singular faces and
indices in every run**, harmonic fields (120 runs) within `2e-13`, the same vertices outside the
cone. `exact` works on the input mesh and has no triangulation.

`tests/representation.rs` (release, jittered grid, jittered icosphere, sheared torus; Delaunay
flipping and Ruppert refinement on both): connectivity equal; stiffness, lumped and Galerkin
mass matrices differ by 0 (bitwise); heat distances from two vertex sources differ by 0, from a
point source by at most 1.8e-15; the Steiner trees have equal length and Steiner points;
adaptive refinement on the icosphere (1 110 vertices, 17 steps) has the same estimate and a
solution difference of 0 at the input vertices.

## Verification

- `tests/exact.rs`: the exact distance equals the Euclidean distance on flat convex domains
  (`1e-10`); bends around a hole at boundary vertices exactly as the shortest paths of the
  plane region (Dijkstra on the visibility graph, an annulus, `1e-10`); on a cube the corners
  are at `1, √2, √5` and points of the faces at their unfoldings (`1e-12`); the distance is that
  of the surface: the same at the input vertices on the input mesh, on the exported iDT and on
  the exported iDR (`1e-9` of the largest distance; icosphere, torus, bumpy grid, sheared
  parallelogram, camelbad); it lies between the straight line in space and the shortest path
  along edges; several sources give the minimum.
- `tests/heat.rs`: close to Euclidean on a flat square (iDT and iDR, mean relative error below
  3 %); a mirror-symmetric grid gives a symmetric distance (`1e-10`); swapping source and target
  on an icosphere changes the value by less than 1 %; on a sheared torus input > 5 × iDT and the
  size-bounded iDR below 2 %; the factorization is reused bit for bit; a vertex given as a point
  of the input mesh gives the same result; values at inserted vertices agree with `value_at` at
  their place; boundary conditions; errors.
- `tests/fem.rs`: the stiffness matrix of a flat intrinsic triangulation (input, iDT, iDR, iODT
  of a jittered grid and a sheared parallelogram) is exact on linear functions at interior
  vertices (`1e-11`) and its Dirichlet energy is `|∇f|²·area`; the Galerkin mass is symmetric,
  positive definite, its rows sum to the lumped mass; the 4-1 subdivision keeps area and angles
  and its lengths are those of the midpoints in space; the Lanczos extreme eigenvalues agree
  with a dense Jacobi eigensolver to `1e-8` (identity and lumped mass) and the Galerkin ones lie
  between the lumped ones and four times them; the pinned Laplacian solves exactly.
- `tests/steiner.rs`: the Fermat point sees the sides under 120° and its tree has the length
  `S² = (a²+b²+c²)/2 + 2√3 A` (`1e-12`); on an icosphere, a cube, a sheared torus and spot's
  control mesh, for the intrinsic and extrinsic MST and Steiner trees: the edges form a
  spanning tree, the cut is one flat disk (Euler characteristic 1, no defect inside, turning
  `2π`, consistent layout), its boundary is twice the tree, every traced path has the length of
  its edge (`1e-8`); on a torus the cut has Euler characteristic −1; the intrinsic trees are not
  longer than the MST of the same triangulation, which is not longer than the extrinsic MST;
  with boundary (a grid, an annulus) the cut is a disk.
- `tests/amr.rs`: on a flat disk bounded by a regular polygon on the unit circle the harmonic
  disk map is a rotation, reproduced exactly, and the estimator vanishes (not on a curved
  disk); the estimator of a smooth solution decreases under uniform refinement; the Green's
  function has mean zero and solves its system, one backward Euler step conserves the
  integral; the marked faces carry the bulk fraction and are the largest; `refine_with` with
  nothing marked is `refine` bit for bit, and every marked face is refined; the adaptive
  solution converges to a fine uniform reference with fewer vertices than uniform refinement
  on a localized problem; the disk map is injective on the iDT (Tutte) and after adaptive
  refinement.
- `tests/vector.rs`: on a flat mesh the transport is the change of frame, every face has
  holonomy zero, and the harmonic vector field is the scalar harmonic function in each
  Cartesian component; the 3D output of a flat field is the field in the plane; both connection
  Laplacians are Hermitian, the finite element one positive semidefinite on every
  triangulation, the finite difference one on the iDT; face holonomies add up to `2πχ` (closed)
  or the interior defects (with boundary); smoothest `n`-fields on spheres have index sum `2n`,
  and on a round sphere the minimal number of singularities; the vector maximum principle holds
  on the iDT; boundary values are interpolated along the boundary (a constant field on a flat
  domain stays constant). `tests/vector_pins.rs`: hashes of every bit the tangent fields
  compute on signposts (9 triangulations), recorded before the code became generic;
  `tests/vector_representation.rs`: the same on integer coordinates (above). Unit tests of
  `vector/knoeppel.rs`: Knöppel's mass entry and
  diagonal, the flat limit (cotan matrix), the nine local entries against a quadrature of the
  model; of `vector/eigen.rs`: the Jacobi eigensolver of the Rayleigh–Ritz step.
- Formula tags (`scripts/check-formulas.sh`): `app-heat-gradient`, `app-heat-divergence`,
  `app-fem-galerkin-mass`, `app-fem-midpoint-subdivision` are machine-checked
  (`lean/proofs/DdgProofs/IntrinsicApps/HeatFem.lean`), as are `app-vector-face-holonomy`,
  `app-vector-fd-connection-laplacian`, `app-amr-map-distortion`
  (`lean/proofs/DdgProofs/Intrinsic/AppsAmrVector.lean`); `app-steiner-fermat-point`,
  `app-vector-transport`, `app-vector-fem-connection-laplacian`, `app-vector-singularity-index`,
  `app-vector-maximum-principle`, `app-amr-residual-estimator`,
  `app-amr-bulk-marking` are cited (`lean/proofs/CITED.md`); the stiffness matrix is
  `cotan-fem-stiffness`, the Laplacian `cotan-laplacian`. As Lebesgue integrals over a triangle
  (`lean/proofs/DdgProofs/TriangleIntegral/`): the Galerkin mass entries, the stiffness entries,
  the area-coordinate rule and slice formula of `knoeppel.rs` (`fem-area-coordinate-rule`) and
  the element norms of `amr::estimate` (`app-amr-element-norms`) and its side fluxes
  (`app-amr-edge-flux`) are machine-checked.

## Reproducing the tables

```text
cargo build --release -p ddg-intrinsic-apps --examples
B=target/release/examples/apps_study
ls MESHES | xargs -P 16 -I{} sh -c "$B steiner {} > runs_tmp/study-steiner-\$(basename {} .obj).log"
ls MESHES | xargs -P 16 -I{} sh -c "$B fem {} > runs_tmp/study-fem-\$(basename {} .obj).log"
for m in MESHES; do $B heat $m > runs_tmp/study-heat-$(basename $m .obj).log; done  # alone: timings
```

Sec. 5.4 and 5.5, one process at a time (the timings are single-thread):

```text
B=target/release/examples
$B/amr_study models/spot/spot/spot_control_mesh.obj green     # also heat, --squares; bob, blub
$B/amr_study jittered_disk disk                                # also bumpy_disk
$B/vector_study smoothest MESHES...
$B/vector_study harmonic sheared:12:4:0.45:2 ... MESHES...
cargo bench -p ddg-intrinsic-apps --bench heat
```

The meshes: `data/models/{spot,bob,blub}/…` (triangulated and control meshes),
`data/optimal-cones/optimal-cones-data/data/{camelbad,camelbadnonuniform,camel_noise1,camelshotnoise,bunny,bump,face}.obj`,
and the generated `icosphere_jitter`, `sheared_torus`, `sheared_square`.
