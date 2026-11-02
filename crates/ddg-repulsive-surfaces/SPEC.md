# Repulsive Surfaces — contract

The contract between the Rust crate (`crates/ddg-repulsive-surfaces`), the Lean oracle
(`lean/oracle`), the Lean proofs (`lean/proofs`), the Python binding (`crates/ddg-py`) and the
Blender add-on for

> Chris Yu, Caleb Brakensiek, Henrik Schumacher, Keenan Crane. *Repulsive Surfaces.* ACM
> Transactions on Graphics 40(6), Article 268, 2021. <https://doi.org/10.1145/3478513.3480521>
> (`../RepulsiveSurfaces/RepulsiveSurfaces.pdf`, "the paper"; equation, section and figure
> numbers are the paper's).

**Reference.** `../RepulsiveSurfaces/repulsive-surfaces` at `c51b464` (MIT, "Copyright (c) 2019
Nick Sharp"; "the reference", cited as `ref: <file>:<line>`). It may be read and ported; ported
code keeps the notice in `LICENSES/repulsive-surfaces-MIT.txt`, and the crate credits all four
authors of the paper. Its submodules (`deps/geometry-central` `00d09a5`, `deps/polyscope`
`feaf66b`, `deps/libgmultigrid` `8d01391`) are empty in the tree; they are fetched only for the
black-box build (WP18), outside this repository, and their licences go into `LICENSES/` before
anything of theirs is ported. A survey of the reference against the paper (2026-11-01) found the
places where the two differ; §11 says in each case what we do and why.

**Sibling crate.** `ddg-repulsive-curves` (the curves paper, same group, MIT) already has the
generic machinery this paper needs: `Bvh<D>`/`Aabb<D>`/`Element<D>`, `BlockClusterTree<N>` with a
`ClusterKernel<N>` trait and the upward/far/downward passes of App. B, `Lbfgs`/`Ncg`,
`ImplicitSurface` with sphere/torus/plane/unions, `Obstacle` (triangle BVH, closest point), the
scalar CCD primitives (`segment_segment_advance`, `point_triangle_distance`, …) and `make_pool`.
All are public. This crate **depends on `ddg-repulsive-curves`** and uses them (D30); it does not
copy them and does not edit that crate. A missing generic piece is written here.

Conventions of `CLAUDE.md` apply: `f64` everywhere, stored indices `u32`, bad input gives an
error enum naming the element and never a panic, nothing is normalized, clamped or repaired
silently, every solve returns a report, the thread count is an option (explicit pool, never the
global `rayon` pool; results bitwise independent of the thread count), formula tags
`[formula: rs-…]` (§12). If an implementation or a test contradicts this file, stop and report
it; do not quietly follow either side.

## 1. The problem

For an immersion `f : M → R³` of a compact surface (several components, boundary and
non-orientable components allowed), the tangent-point energy (Eq. 2, 3, in the reference's
two-exponent form) is

```
ℰ(f) = ∬_{M²} |N(x)·(f(x) − f(y))|^α / |f(x) − f(y)|^β  dx_f dy_f
```

with `N` the unit normal (the sign is irrelevant). The paper uses `α = p`, `β = 2p`, `p = 6`.
The differential has order `2s`, `s = (β − 2)/α` (`= 2 − 2/p` for `β = 2p`), `σ := s − 1`. The
energy is a barrier against self-contact when `β − α > 4` (`p > 4`); for `β − α ≤ 4`
("subcritical", Fig. 26) it is finite on intersecting surfaces but still discourages overlap.

We minimize `F = ℰ̂ + Σ penalties` over the vertex positions subject to constraints
`Φ(f) = 0`, by projected gradient descent in the discrete fractional Sobolev metric
`A = B + B⁰` (Sec. 2.4, 3.2), with Armijo line search, optional continuous collision detection,
constraint projection, and dynamic remeshing after every step (Sec. 6). Acceleration: Barnes–Hut
for the energy and its gradient (Sec. 4), hierarchical matrices for products with `A` (Sec. 5.1,
App. A, B), GMRES preconditioned by `(−Δ)⁻¹ L_{2−s} (−Δ)⁻¹` (Sec. 5.2), Schur complements for
nonlinear constraints and saddle rows for positional ones (Sec. 5.3, 7.2). For the evaluation of
Sec. 8 the same machinery also runs L², H¹, H² projected gradient, AQP, H¹ L-BFGS and BCQN.

In scope (user, 2026-11-01): everything above; every constraint and penalty of Sec. 7; all
applications of Sec. 9 (§14.4); all six comparison methods (§14.2); a CCD guard in the line search
(on by default only if measured cheap, D20); Blender as for curves (§9.3).

## 2. Discrete setting and data model

### 2.1 Surface mesh (D1)

`M = (V, E, F)`, positions `f_v ∈ R³`, faces `S = (i, j, k)` in the caller's vertex order.
`Surface::new(n_vertices, faces)` validates and stores: edges, edge–face incidence, one-ring
fans, connected components (face- and vertex-wise), boundary loops, and the **orientation class**
of each component: `Consistent` (every interior edge is traversed in opposite directions by its
two faces), `Orientable` (consistent after flipping some faces; the caller's order is kept and
the flip set reported) or `NonOrientable` (Klein bottle, Fig. 26).

* Manifold is required: every edge has one or two faces, every vertex star is one fan (disk or
  half-disk). Errors name the element, in the style of `ddg_core::MeshError`
  (`"edge (12, 57) is shared by 3 faces: 4, 9, 31"`).
* Orientation is **not** required: `ddg_core::TriMesh` rejects non-orientable meshes, so this crate
  has its own validated type. Nothing in the energy, the metric, the barycenter, pins, area,
  remeshing or CCD depends on the orientation (the kernel uses `|N·v|`, `D_f` is the intrinsic
  P1 gradient, §4.1). Signed volume needs `Consistent` components: a volume constraint or
  penalty on any other component is an error naming an edge whose two faces disagree.
* Triangles only. Blender passes its loop triangles; polygons are never triangulated here.
* Degenerate faces (zero area) are errors at validation of positions; non-finite positions too.
* Unreferenced vertices are an error (they have no area and break the metric).
* Components are numbered by their lowest face index (WP1, oracle). `Orientable { flips }` keeps
  the component's lowest face fixed; `Component::inconsistent_edge()` names the edge reported by
  `VolumeNeedsOrientation`.

Per face: area `a_S = ½|n_S|`, `n_S = (f_j − f_i) × (f_k − f_i)`, unit normal `N_S = n_S/|n_S|`,
barycenter `X_S = (f_i + f_j + f_k)/3`. Edge vectors `e_jk = f_k − f_j`, `e_ki = f_i − f_k`,
`e_ij = f_j − f_i` (not normalized; D2). Vertex (barycentric dual) area `a_v = ⅓ Σ_{S∋v} a_S`.
Cotan stiffness `W` (PSD, `W_ij = −½ Σ cot`), lumped mass `M = diag(a_v)` — computed here from the
face list (they are orientation-free); tests compare with `ddg_core::operators` on orientable
meshes.

### 2.2 Exponents (D3)

`Exponents::new(α, β)` requires finite `α > 1`, and `s = (β − 2)/α ∈ (1, 2)`, i.e.
`α + 2 < β < 2α + 2`. Default `(6, 12)`. `Exponents::subcritical()` is `β − α ≤ 4`; it is allowed
(Fig. 26) and reported in every `EnergyReport`. Derived: `s`, `σ = s − 1 ∈ (0, 1)`, high/low
kernel exponent `2s = 2σ + 2`, preconditioner kernel exponent `6 − 2s = 2(2 − s) + 2`.

### 2.3 Input validation

`validate_positions` checks count, finiteness, zero-area faces. Coincident barycenters of two
distinct faces give `ℰ = ∞`: reported as `NonFiniteEnergy { faces: (S, T) }` at setup; inside a
line search a non-finite trial energy is a rejected trial, never an error.

## 3. The discrete energy

### 3.1 Kernel and energy (Eq. 9, 10, 16)

```
K(X, N; Y) = |N·(X − Y)|^α |X − Y|^{−β}
ℰ̂(f) = Σ_{S ∈ F} Σ_{T ∈ F, T ≠ S} K(X_S, N_S; X_T) a_S a_T
```

over ordered pairs, both orientations, `S = T` excluded, neighbours **not** excluded (paper and
`ref: tpe_all_pairs.cpp:52,65` agree). No `2^p`, no ½. Integer `α`, `β` use integer powers
(bit-reproducible, faster); others `powf`.

### 3.2 Differential

Chain rule through per-face quantities. For face `S = (i, j, k)` and its vertex `i`:

```
∂X_S/∂f_i = I/3
∂a_S/∂f_i = ½ N_S × e_jk                                   (rs-face-area-gradient)
∂N_S/∂f_i = (1/(2a_S)) (e_jk × N_S) N_Sᵀ                  (rs-face-normal-variation)
```

(cyclically for `j`, `k`). With `v = X − Y`, `ν = N·v`, `r = |v|`:
`∂K/∂v = α|ν|^{α−2}ν r^{−β} N − β|ν|^α r^{−β−2} v`, `∂K/∂N = α|ν|^{α−2}ν r^{−β} v`,
`∂K/∂Y = −∂K/∂v` (`rs-kernel-gradient`; `|ν|^{α−2}ν` is continuous at `ν = 0` for `α > 1`).
Whether the code accumulates in `(a, X, N)` or in the reference's `(a, aX, aN)` variables is free;
the result is the exact gradient, checked by finite differences (§13.2).

### 3.3 Invariances and scaling

`ℰ̂` is invariant under rigid motions and under flipping any face's orientation;
`ℰ̂(cf) = c^{4 + α − β} ℰ̂(f)` (`rs-energy-scaling`). Hence `Σ_v g_v = 0`, `Σ_v f_v × g_v = 0`,
`Σ_v f_v·g_v = (4 + α − β) ℰ̂` (`rs-gradient-identities`), tested to round-off.

### 3.4 Barnes–Hut energy (Eq. 17; D4, D5)

A BVH over faces (`Bvh<3>`, element point `X_S`, bounds the triangle's box, mass `a_S`; leaf
size 8). Node `I` stores `A_I = Σ a_S`, `X_I = Σ a_S X_S / A_I`, the AABB of its triangles'
**vertices**, the radius `r_I` = distance from `X_I` to the farthest box corner, and (for the
metric) `P̄_I = Σ a_S N_S N_Sᵀ / A_I`. The traversal is per **leaf** `ℓ` (`ref:
tpe_barnes_hut_0.cpp:100`, D5): node `C` is admissible for leaf `ℓ` iff
`max(r_ℓ, r_C) < θ · d_box(ℓ, C)` with `d_box` the distance between the two boxes (zero when they
overlap, so overlapping boxes are never admissible). Admissible: every face `S ∈ ℓ` adds
`K(X_S, N_S; X_C) a_S A_C`. Inadmissible leaf × leaf: exact face sum, `S ≠ T`. Default
`θ = 0.5`; `θ = 0` is the exact energy (every node inadmissible). The admissibility lists are
built once per step and frozen during the line search (paper Sec. 4.1); trials refresh only the
face data and the node aggregates `A`, `X`, `P̄` (`refit`); boxes and radii stay those of the
build. Leaves may be lumped too (as the reference). Measured (WP3, 10k-face torus): relative
energy error 8.4e-2 / 2.3e-2 / 4.1e-3 / 2.1e-4 at `θ` = 1 / 0.5 / 0.25 / 0.125 (at least
`O(θ²)`); 40k faces, energy + gradient 0.57 s on 1 thread, 0.15 s on 8. The default `θ` is
revisited in §14.3 (2.3 % energy error at 0.5 is large for the line search's comparisons).

**Gradient (D4).** The gradient is the **exact gradient of `Ẽ` for the frozen lists**, including
the terms through the cluster aggregates `A_C`, `A_C X_C` (pushed down the tree to the member
faces; `ref: optimized_cluster_tree.cpp:1246–1311`). The paper's Sec. 4.2 formula keeps only the
derivatives through `S` and would not be the gradient of the energy the line search compares.
`rs-bh-aggregate-adjoint`: for `A_C = Σ a_T`, `M_C = Σ a_T X_T`, `∂Ẽ/∂a_T = Σ_{C ∋ T} ∂Ẽ/∂A_C`
and `∂Ẽ/∂(a_T X_T) = Σ_{C∋T} ∂Ẽ/∂M_C`.

### 3.5 Other energies

`energy()` (§9.1) evaluates the tangent-point energy (exact or Barnes–Hut) plus penalties; the
session (WP5) dispatches `Energy::Coulomb` and `Energy::Willmore` to WP7's code.
`Energy::TangentPoint` (default). `Energy::Coulomb { alpha }` (Fig. 2, evaluation only, exact
dense, small meshes): `Σ_{v ≠ w} a_v a_w |f_v − f_w|^{−α}` over ordered vertex pairs with dual
areas (the paper's "penalizing the distance between all pairs of vertices"; the reference's
variant has no areas, a fixed exponent 1, half the gradient and erases the other gradients, D24).
`Energy::Willmore` replaces `ℰ̂` by §5.3's Willmore term (the reference's `method willmore`, used
with the H² metric for Fig. 4).

## 4. Inner product

### 4.1 Operators on vertex functions

`ū(S) = (u_i + u_j + u_k)/3` (`U u = a ⊙ ū`, the paper's `U`). The P1 gradient (D2)

```
D_f u(S) = (1/(2a_S)) N_S × (u_i e_jk + u_j e_ki + u_k e_ij)          (rs-df-gradient)
```

with **un-normalized** edge vectors: it is the unique `g ⊥ N_S` with `g·(f_b − f_a) = u_b − u_a`
on the triangle. The paper's "unit edge vectors" (Sec. 3.2.1) is a misprint (the unit-edge form
is off by 0.85 relative on a random triangle; `ref: hs_operators.h:104–115` uses the plain
edges). `V = diag(a) D_f` (3 rows per face). `D_f 1 = 0`, `Ū 1 = 1`.

### 4.2 The three operators (Eq. 11–13), over ordered pairs `S ≠ T`, `ρ = |X_S − X_T|`

```
uᵀ L_τ v = Σ_{S≠T} (ū_S − ū_T)(v̄_S − v̄_T) a_S a_T ρ^{−(2τ+2)}                    (fractional, order 2τ)
uᵀ B   v = Σ_{S≠T} ⟨D_f u(S) − D_f u(T), D_f v(S) − D_f v(T)⟩ a_S a_T ρ^{−2s}
uᵀ B⁰  v = Σ_{S≠T} (ū_S − ū_T)(v̄_S − v̄_T) K₂(S, T) a_S a_T ρ^{−2s},   K₂(S,T) = |N_S·(X_S − X_T)|² ρ^{−4}
```

`B⁰` uses the fixed `K₂` (α = 2, β = 4) whatever the energy's exponents (paper, reference). As an
ordered sum it equals the sum with the symmetric kernel `½(K₂(S,T) + K₂(T,S))`
(`rs-low-order-symmetrisation`). The dense assembly loops over face pairs and the at most six
vertices of the pair (Alg. 1; the paper's "distinct pairs" there means ordered pairs, consistent
with Eq. 11 and App. A).

### 4.3 Kernel-matrix form (App. A)

For a symmetric kernel matrix `H` with zero diagonal and positive `a`:

```
Σ_{S,T} (ū_S − ū_T)(v̄_S − v̄_T) a_S H_ST a_T = 2 (Uu)ᵀ (diag(a)⁻¹ diag(H a) − H) (Uv)     (rs-kernel-matrix-form)
```

and the same with `V` and vector-valued `D_f u` (sum over the three components). The
hierarchical product computes `H a` with **the same compressed `H`** as the product (`ref:
optimized_bct.cpp:936–1019`), so `L·1 = 0` holds exactly for the compressed operator too
(`rs-hierarchical-null-space`); `B·1 = 0` holds by `D_f 1 = 0`.

### 4.4 The metric `A` (D6, D7)

```
A = w_TP (B + c_lo B⁰) + w_W (W M⁻¹ W) + O_diag,        A₃ = A ⊗ I₃
```

* `w_TP` the tangent-point weight (default 1; the reference scales its metric by it, `ref:
  optimized_bct.cpp:332`). `c_lo ∈ {0, 1}`: `low_order: false` drops `B⁰` (Sec. 9.2.4, Fig. 26;
  the reference's `disable_near_field`, which drops the whole low-order term, not only its near
  field — the energy is unchanged).
* `w_W W M⁻¹ W` only with a Willmore penalty of weight `w_W` (Sec. 7.3.5). The reference adds it
  unweighted and only inside GMRES; we weight it (it is the Hessian-like term of `w_W E_W`) and,
  being sparse, also add it to the preconditioner's Laplacian part (D7: measured; reference mode
  available).
* `O_diag` the mesh-obstacle term (§5.3, D23): for each obstacle `𝒪` with weight `w_𝒪`, the
  diagonal `2 w_𝒪 Uᵀ diag(Σ_{T∈𝒪} h(S,T) a_T / a_S) U` and its `V` analogue — the metric of
  mesh ∪ obstacle with the obstacle fixed. Every obstacle contributes (the reference keeps only the
  last one added, `ref: surface_flow.cpp:47`).
* No mass term (paper, reference). Kernel of `A₃` (WP2): with `B⁰` on, only the **global**
  translations (`B⁰` couples components: `ū_S − ū_T ≠ 0` across components for per-component
  constants); with `low_order: false`, the translations of each component. The constraints remove
  it (§6.1, D8).

### 4.5 The comparison metrics (Sec. 8.2.1, D25)

`L2`: `M ⊗ I₃` (paper: "the mass matrix"; the reference uses the identity — option
`L2Kind::Identity`). `H1`: `W ⊗ I₃`. `H2`: `W M⁻¹ W ⊗ I₃`. Each solved directly (sparse LDLT) with
all constraints in the saddle (as the reference).

## 5. Constraints and penalties (Sec. 7)

### 5.1 Constraints

Two kinds (Sec. 7.2): **positional** rows live inside the operator `A′` (and inside the
preconditioner's Laplacian), **nonlinear** ones go through the Schur complement.

| constraint | `Φ` | Jacobian row(s) | kind | projection |
|---|---|---|---|---|
| `Barycenter { component, target }` | `Σ_v a_v f_v / Σ_v a_v − X₀` over a component | `a_v / A_c ⊗ I₃` (areas frozen, D8) | positional, 3 rows | exact translation of the component |
| `Pins { vertices, targets }` | `f_v − f_v⁰` | `e_v ⊗ I₃` | positional, 3 rows each | exact reset |
| `Area { component, target }` | `Σ a_S − A₀` | `Σ_{S∋v} ½ N_S × e_opp` | Schur, 1 row | Newton (§6.2) |
| `Volume { component, target }` | `⅙ Σ f_i·(f_j × f_k) − V₀` | `⅙ Σ_{S∋v} f_j × f_k` (= ⅓ Σ a_S N_S on a closed fan) | Schur, 1 row | Newton |

* `component: None` means the whole mesh for area/volume. For the barycenter it follows
  `BarycenterMode` (D8, revised after WP2): `Auto` (default) uses **one global barycenter** when
  `B⁰` is on — the kernel of `A₃` is then the global translations only, and components move
  relative to each other through the metric itself (the true `Hs` gradient; also what the
  reference does) — and **per-component barycenters plus the per-component L² mean-motion shift**
  (paper Fig. 8, Sec. 7.1.1) when `low_order: false`, where each component's translations are in
  the kernel. `Global` and `PerComponent` force one or the other (`PerComponent` without the shift
  freezes the components' barycenters: a real constraint, e.g. for packing scenes). In `Global`
  mode on several components the preconditioner's `W` keeps per-component constants in its kernel:
  `L′` then uses `ε_multi = 1e-4 · tr W / tr M` (WP6: the reference's `1e-2`, made scale-invariant,
  costs up to 6× more GMRES iterations; at `1e-4` Global costs about what PerComponent does). WP6 measures `Auto` against `PerComponent` + shift on Hopf link, Borromean
  rings and the keys/double-torus-and-ring scene (GMRES iterations, energy per step). **Measured
  (WP6):** Global and PerComponent directions differ by 46 % (Hopf) and 67 % (Borromean): freezing
  the components' barycenters is a real constraint. The paper's L² mean-motion shift is not
  scale-consistent: on the Borromean rings it made the direction a rigid translation (the shift
  was 2e4 × ‖x‖). Hence `Auto` = `Global` with `B⁰` on (default), and with `low_order: false`
  `PerComponent` **without** the shift; `barycenter_shift: true` stays available as the paper's
  option (Fig. 8, Fig. 18 free barycenter), documented as unit-inconsistent.
* A barycenter together with pins on a component is not linearly dependent (only a fully pinned
  component is), but it constrains more than the pins alone; the projection iterates
  translate-then-reset until both hold (WP2). The paper: "if any pins are used, a barycenter
  constraint is unneeded".
* `target: None` = the initial value. A barycenter target with `component: None` applies to every
  component (e.g. concentric envelopes); use one constraint per component for different targets.
  Volume on an open component is the signed volume of the
  cone over the origin (origin-dependent; used by the eversion of Fig. 20 with a fixed barycenter).
* **Volume Jacobian (D9).** We use the exact gradient. The reference's row is `Σ a_S N_S = 3∇V`
  (`ref: total_volume.cpp:30`): the descent direction is unchanged but its one Newton step removes
  only a third of the volume error per step. `rs-volume-gradient` proves the factor.
* The translation null space of `A₃` must be removed on every component (barycenter or ≥ 1 pin);
  otherwise `TranslationNotFixed { component }`. Rows that are linearly dependent (e.g. a
  barycenter and pins on a component, three pins that fix it already) are detected by the rank
  check of the curves crate's pattern and the redundant ones dropped with a report (D10).
* `BarycenterShift` (Sec. 7.1.1 "Barycenter Motions"): with the barycenter constraint used only to
  fix the null space, after the solve add per component the `a_v`-weighted mean of the L²
  differential to the direction (as printed: it mixes an L² differential with an `Hs` direction;
  the reference does it once for the whole mesh and then skips all positional projections,
  including pins — we skip only the barycenter's).
* **Schedules (D11).** `Schedule { constraint, factor, offset, steps }`: the target moves linearly
  from `T₀` to `factor·T₀ + offset` over `steps` steps (the reference's `constrain area 3 200 1.5`).
  Moving pins: `PinSchedule { displacement, steps }` (reference `constrain vertices … move`). The
  reference raises the target after the step's projection; we raise it before (the step projects
  onto this step's target) — reference mode available.

### 5.2 Penalties (D24: exact gradients of the printed values)

Each penalty is `w · P(f)` with a value and an exact gradient; finite differences test every one.

| penalty | `P(f)` | source / note |
|---|---|---|
| `AreaDeviation` | `(A/A₀ − 1)²` | Sec. 7.3.1. Reference: value same, gradient has a spurious `a_v` factor |
| `VolumeDeviation` | `(V/V₀ − 1)²` | Sec. 7.3.1. Reference gradient 3× |
| `TotalArea { power }` | `A^power` (1: Fig. 9 minimal surfaces; 2: the reference's `minimize area`) | reference gradient ½× for power 2 |
| `TotalVolume { power }` | `V^power` (signed) | reference gradient 3/2× |
| `SquaredError` | `Σ_v |f_v − f_v⁰|²` (Fig. 21–22) | targets carried through remeshing (D17); the reference's gradient is wrong with > 1 thread |
| `BoundaryLength { target }` | `(L − Σ_{e∈∂M} ℓ_e)²` | Sec. 7.3.4. Reference returns `L_∂ − L` unweighted, unsquared |
| `BoundaryCurvature { dual_length }` | `Σ_{v∈∂M} θ_v² / ℓ_v` (`dual_length: false`: no `ℓ_v`, the reference) | Sec. 7.3.4; `θ_v` the turning angle, `ℓ_v` half the two boundary edges |
| `Willmore` | `fᵀ W M⁻¹ W f` (Sec. 7.3.5; interior rows; boundary rows of `W f` dropped as the reference) | gradient includes the variation of `W` and `M` (derived, FD-tested); the reference's value is `¼` of this |
| `ImplicitObstacle { surface, power }` | `Σ_v q_v d(f_v)^{−power}` | Sec. 7.3.3; `d ≤ 0` gives `+∞` (rejected trial) |
| `ImplicitAttractor { surface, power, vertices }` | `Σ_v q_v |d(f_v)|^{power}` | Sec. 7.3.3; the reference picks attracted vertices by UV `u > 1e-10`, we take an explicit vertex set (`None` = all) |
| `MeshObstacle { obstacle, kind }` | §5.3 | Sec. 7.3.2 |
| `PointCloudObstacle { points }` | `Σ_S Σ_p K(X_S, N_S; y_p) a_S` (one-sided TP, unit point weights) | Fig. 23; Barnes–Hut over a point BVH |

`q_v` is the quadrature of the implicit potentials (D21): `Quadrature::DualArea` (`q_v = a_v`, the
default — a discretization of `∫ d^{−p} dA`, refinement-stable) or `Quadrature::PerVertex`
(`q_v = 1`, the reference; scene import uses it so the reference's tuned weights carry over).
Implicit surfaces come from the curves crate (`Sphere`, `Torus`, `Plane`, unions, …); this crate
adds `Cylinder { point, axis, radius }` (infinite) and the reference's `y`-axis torus as a
`Torus` with a frame.

### 5.3 Mesh obstacles (D23)

The paper prints a Coulomb-like potential `E_obs(x) = Σ_{S∈F_𝒪} |f_𝒪(S) − x|^{−p} a_𝒪(S)`
(Sec. 7.3.2); the reference implements a **two-sided tangent-point cross energy**
`Σ_{S∈M} Σ_{T∈𝒪} K(X_S, N_S; X_T) a_S a_T + Σ_{T∈𝒪} Σ_{S∈M} K(X_T, N_T; X_S) a_T a_S`
with Barnes–Hut (same `θ`) plus the metric term `O_diag` (§4.4). We implement both:
`ObstacleKind::TangentPoint` (default; it is what produced the walnut, lung and hamburger scenes,
and it keeps the barrier property) and `ObstacleKind::Coulomb { power }` (the printed formula,
integrated over the mesh with `a_S` at barycenters: `Σ_{S∈M} Σ_{T∈𝒪} a_S a_T |X_S − X_T|^{−p}`, mesh barycenters against obstacle-triangle centroids; no metric term). `penalties::evaluate` takes the energy's `Exponents` and `EnergyEval` (the TP obstacle uses the same exponents and `θ`; `θ = 0` for `Exact`). In the line search the session keeps an `obstacle::CrossEvaluator` per obstacle (built per step, refit per trial, lists frozen), as for the energy (WP7). The obstacle code carries the existing ids (`rs-kernel`, `rs-energy-gradient`, `rs-bh-*`); no separate ids. Obstacles are static
triangle soups (`Obstacle` of the curves crate), need not be manifold. A point-cloud obstacle (one-sided, unit weights) also has the metric term `O_diag`
in the reference (zero obstacle projectors, so `K₂(T,S) = 0`); we add it in `reference_compat` only (`CrossEvaluator::metric_diagonal`;
the preconditioner lacks it and default-mode GMRES then needs ≈ 1600 iterations per step on the Reconstruction scene).

## 6. Descent direction and projection

### 6.1 Constrained direction (Eq. 14, 18–20; Sec. 7.2)

`g = dF(f)` (energy + penalties, `3|V|`). With the positional rows `C_Q` inside
`A′ = [[A₃, C_Qᵀ], [C_Q, 0]]` and the nonlinear rows `C_N` (`k` of them), the direction is

```
x = A′⁻¹ g + A′⁻¹ C_Nᵀ S⁻¹ C_N A′⁻¹ g,     S = −C_N A′⁻¹ C_Nᵀ   (k × k, dense)        (rs-schur-direction)
```

(restricted to the position block; `C_N` padded with zeros on the multiplier columns). `k + 1`
solves with `A′`; `A′⁻¹ C_Nᵀ` is kept for the projection. `C x = 0` up to the solver tolerance
(reported). The step is `f ← f − τ x`.

### 6.2 Constraint projection (Eq. 15, 21; D12)

After the line search: positional constraints exactly (translate each constrained component, reset
pins). Nonlinear ones by the Newton step `h = −A′⁻¹ C_Nᵀ S⁻¹ Φ(f)` (position block; with `S` as in
§6.1 this gives `C_N h = Φ`, `rs-schur-correction`), applied as `f ← f − h`, reusing `S` and
`A′⁻¹C_Nᵀ` of the step (old Jacobian; the paper: "repeated several times if the violation is not
sufficiently small"), each round followed by the positional projection (Newton moves the
barycenter). Default: repeat while `|Φ_r| > tol_r` (`1e-8 |target|`, or `1e-8` for a zero target)
up to 5 times. The reference does exactly one step; option `newton_steps: Fixed(1)`.

**Safeguard (found by the oracle's selftest).** The frozen-Jacobian iteration can diverge: on a
near-sphere `∇A ≈ (2/r)∇V`, so `S` is nearly singular when both area and volume are constrained,
and one round moved vertices by a mean edge length (violations of 1e80 after a few rounds). A
round is kept only if it lowers `max_r |Φ_r|/tol_r`; otherwise it is undone, the rounds stop and
the step reports `ProjectionStatus::NotImproving`. Before the solve, `S` is equilibrated and
rows whose pivot falls below `1e-10` relative are dropped for this step (reported in
`dropped_rows`, D10), so area and volume together on a sphere keep the better-conditioned row
instead of producing a huge correction. A violation above `tol` after the cap is reported, not an
error.

### 6.3 Solver for `A′` (Sec. 5, D13, D14)

* **Dense** (`Backend::Dense`, small meshes and the oracle): assemble `A₃` (Alg. 1), saddle with
  `C_Q`, LU (or Cholesky of the reduced system), direct Schur.
* **Iterative** (`Backend::Hierarchical`, the paper): GMRES on `A′` with the products of §8.3 and
  the preconditioner

  ```
  P⁻¹ r = L′⁻¹ · (blockdiag(L_{2−s} ⊗ I₃, 0)) · L′⁻¹ r,    L′ = [[(W + ε M) ⊗ I₃, C_Qᵀ], [C_Q, 0]]   (rs-fractional-preconditioner)
  ```

  (the multiplier part is zeroed between the two solves; `L′` factored once per step by sparse
  LDLT).

  **The fractional factor (D32, after WP2).** The paper/reference evaluate `L_{2−s}` as the face
  form of §4.2 (`ū_S − ū_T`, face pairs). Face averaging nearly annihilates high-frequency vertex
  functions, so this `L_{2−s}` is a poor stand-in for `M Δ^{2−s}` there: condition numbers
  12.6 / 68 / 224 against `M Δ^{2−s}` on 80 / 320 / 1280 faces, and GMRES to `1e-4` needs 38 / 74 /
  > 200 right-preconditioned iterations on icosphere 320 / 1280 / torus 1500 and stalls at `4e-2` on a
  bumpy disk. Left-preconditioned GMRES (the reference) stops after 13–19 iterations but its true
  residual is then 0.1–0.2: the reference's solves are much less accurate than its tolerance
  says. The **vertex-pair** form
  `uᵀ L^V_{2−s} u = Σ_{v≠w} (u_v − u_w)² a_v a_w |f_v − f_w|^{−(2(2−s)+2)}` (dual areas, ordered
  pairs) has condition 1.3–2.1 against `M Δ^{2−s}` and GMRES needs 5 / 5 / 7 / 27 iterations on the
  same meshes. Default: `FractionalForm::VertexPairs` (hierarchical on a vertex BVH, WP6);
  `FractionalForm::FacePairs` is the reference mode. Confirmed by WP6 at 2k–100k faces: vertex
  pairs, right GMRES to `1e-4` in 4–16 iterations (Borromean rings 37); face pairs, right GMRES does
  not converge in 200 on most meshes, left stops after 4–16 with true residuals 0.2–2.2.

  **The augmented form (D33, after WP6).** As printed, `P⁻¹ = L′⁻¹ blockdiag(L_{2−s}⊗I₃, 0) L′⁻¹`
  is singular (the middle factor zeroes the multiplier block), so its range cannot represent the
  multiplier part of the solution: whenever pins or per-component barycenters carry a force, right
  GMRES stalls (Hopf link 8.6e-4 after 200 iterations; 3 pins on a torus 0.12) and left GMRES stops
  with true residuals 0.04–0.11. Default:

  ```
  P⁻¹ = L′⁻¹ · [[L_{2−s}⊗I₃, C_Qᵀ], [C_Q, 0]] · L′⁻¹  =  L′⁻¹ + L′⁻¹ blockdiag((L_{2−s} − (W+εM))⊗I₃, 0) L′⁻¹
  ```

  (same cost: two back-substitutions and one fractional product); symmetric, nonsingular, and
  `C_Q`-consistent: `C_Q (P⁻¹[r; b])_x = b`. The printed form stays as `PreconditionerForm::Printed`
  (reference mode). `rs-fractional-preconditioner` names this form (proved, `RS.aug_precond_props`). Not added to `L′`: the Willmore term (D7: with it GMRES needed 154 iterations or diverged,
  5–20 without) and the obstacle diagonal (the reference adds its `fr` part; not implemented, §16). `ε = ε_rel · tr(W)/tr(M)` (scale-invariant; default `ε_rel = 1e-10`; the reference uses
  absolute `1e-8`, and `1e-2` with several components because its global barycenter row leaves
  per-component constants in the kernel — with per-component barycenters we need no large `ε`).
  The order bookkeeping is `rs-preconditioner-order` (`(−Δ)⁻¹ (−Δ)^{2−s} (−Δ)⁻¹ = (−Δ)^{−s}`, with
  stiffness matrices the mass factors cancel — proved for the matrix identity
  `W⁻¹ (M Δ^{2−s}) W⁻¹ = Δ^{−s} M⁻¹` with `Δ = M⁻¹W`); the approximation of `B` by it is cited.
* **GMRES** (D13): restarted, **right** preconditioned (it minimizes the true residual; the
  reference's Eigen GMRES is left preconditioned and stops on the preconditioned residual —
  `GmresSide::Left` reproduces it), initial guess 0, restart 30, relative true-residual tolerance
  `1e-4`, at most 200 iterations by default (the reference: `2n`), modified Gram–Schmidt, all
  `x`/`y`/`z` columns as one `|V| × 3` block where the operator allows (App. B). Not converging is
  a reported condition (`LinearReport`), not an error, unless the residual is above `1e-1`
  (`SolverDiverged`). **Recovery:** a diverged solve of the hierarchical back end is repeated by right
  GMRES with the augmented preconditioner on vertex pairs (restart 60, ≥ 800 iterations, rebuilt preconditioner),
  then by a dense LU of the metric (≤ 6000 vertices); the escalation is listed in `LinearReport::fallbacks`
  and `SolverDiverged` is returned only when both fail. With `Backend::Auto` the dense path is used up to `dense_max_faces`
  (default 1000: WP6 measured the hierarchical path faster from about 500 faces).

## 7. Time stepping

### 7.1 One step

1. Build the face BVH, the Barnes–Hut lists and the block cluster tree on the current positions.
2. `g`; with `BarycenterShift`, its per-component weighted means.
3. Factor `L′`; the direction `x` (§6.1); the shift.
4. Line search along `−x` (§7.2), optionally capped by CCD (§7.3).
5. Projection (§6.2). Schedules advance.
6. Remeshing (§10) unless off; the per-vertex data (pins, targets, attractor sets, squared-error
   targets) follow the vertex map.
7. Stop rules (§7.4).

`Session` keeps the mesh, positions and per-vertex data between steps (Blender live mode), is
cancellable between steps and between line-search trials.

### 7.2 Line search (D15)

Backtracking on `F` along `f − τ x`, halving.

* `ArmijoRule::Standard { c = 1e-2 }`: accept if `F(f − τx) ≤ F(f) − c τ ⟨g, x⟩` (default).
* `ArmijoRule::Reference { c = 1e-2 }`: accept if `F(f) − F(f − τx) ≥ c τ ‖x‖ cos∠(g, x)`
  (`ref: line_search.cpp:86,103`: the `‖g‖` factor is missing, so the test is not scale-invariant;
  `rs-armijo-reference` proves it equals the standard test with `c' = c/‖g‖`).
* Start (`StartRule`): `Reference` (`τ₀ = 2/‖x‖_F`, absolute units, `ref:
  surface_flow.cpp:60–66`), `MaxDisplacement { lengths }` (`τ₀` moves the farthest vertex by
  `lengths · ℓ̄₀`, scale-invariant), `GrowPrevious { factor, first }`. Provisional default:
  `GrowPrevious { factor: 2, first: MaxDisplacement { lengths: 1 } }`; the default is fixed by the
  measurement of §14.3.
* Stop at `τ < τ_min = 1e-12 · ℓ̄₀/‖x‖_∞` (the reference: `1e-20` absolute): the step is `Stuck`,
  positions unchanged; projection, schedules and remeshing still run (as the reference), and
  `stuck_limit` (default 5) consecutive stuck steps end the run with `Status::Stuck`.
* The Barnes–Hut lists and the BCT are frozen during the search; every trial refits face and
  node data; a non-finite trial value rejects the trial.

### 7.3 Continuous collision detection (D20, not in the paper)

`Ccd { gap }`: before the line search, the largest `τ_max ≤ τ₀` such that the linear motion
`f − t x`, `t ∈ [0, τ_max]`, keeps every vertex–triangle and edge–edge pair of non-adjacent
elements (and mesh–obstacle pairs) at distance ≥ `gap · ℓ̄₀` (default `gap = 0.05`). The test is
conservative advancement (the curves crate's `rc-ccd-advancement`, extended to moving triangles:
relative motion bound `max |u_p| + max |u_tri|`), with a BVH over swept boxes as broad phase. The start becomes
`min(τ₀, ⅔ τ_max)` (as the curves crate). **As built (WP8):** the gap is `gap · ℓ̄₀` absolute; a pair already
closer than the gap at `t = 0` may close to `close_factor · d₀` (default 0.5) and a pair at distance
0 is unconstrained, both counted in the report (`close_pairs`, `touching_pairs`), so the flow never
freezes; the motion bound is the smaller of `max |u_p| + max |u_tri|` and the same bound in the
pair's mean-velocity frame; the broad phase is a parallel Morton box tree (the curves crate's
`Bvh<3>` builds sequentially and has no per-element box test). The collapse guard also refuses a
collapse next to a pair already closer than the gap. Cost (WP8, 21 cases, 5k–50k faces, machine
load 10–18): median CCD time / Barnes–Hut energy+gradient time 0.29 (1 thread), 0.19 (8 threads);
the share of a whole step is measured by WP5. Remeshing can create intersections too (a collapse
moves vertices): the remesher checks its own moves (§10). Default: **off until measured**; WP8
measures the per-step cost on the Fig. 13 set and the number of runs where the energy alone lets
surfaces cross (coarse meshes, Fig. 27); if the median cost is below 10 % of a step, it becomes
the default (the user's decision: "default by measurement"). Recorded in §17.

### 7.4 Stop rules and status

`max_iterations` (default 1000), `time_limit`, `StopRule::{EnergyStall { relative, window },
SobolevNorm { eps }}` (`⟨g, x⟩` is the squared `Hs` norm of the gradient; relative to its first
value). Default `EnergyStall { relative: 1e-7, window: 10 }`. `Status::{Running, Converged,
MaxIterations, TimeLimit, Stuck, Cancelled, Error(e)}`; `run` never returns an error after the
first step (it ends with `Status::Error`), as the curves crate.

## 8. Acceleration (Sec. 4, 5, App. B)

### 8.1 BVH

`Bvh<3>` from the curves crate over faces (§3.4); its split rule (axis cycling, diameter
criterion) differs from the reference's longest-axis midpoint split — the tree shape changes the
approximation only within the admissibility bound (D5). Rebuilt each step (remeshing changes the
faces); refit during the line search. Per-node payload via `Bvh::fold_up`.

### 8.2 Barnes–Hut energy and gradient

§3.4. Deterministic: parallel over leaves, per-leaf partial sums combined in a fixed order.

### 8.3 Block cluster tree and products (Sec. 5.1, App. A, B)

`BlockClusterTree` of the curves crate on the face BVH; `ClusterKernel` implements the paper's
separation `max(r_I, r_J) ≤ χ · d_box(I, J)` (`≤`, `χ = θ = 0.5` by default as the reference), far
kernels from node data (`X_I`, `P̄_I`, `A_I`), near kernels from face data. Three kernel families
on one tree (App. B): `h_hi = ρ^{−2s}` (for `B`, channels `V u`, 3 per input column), `h_lo =
½(vᵀP̄_I v + vᵀP̄_J v) ρ^{−4} ρ^{−2s}` (for `B⁰`; `vᵀP̄v` is the area average of `|P_S v|²`,
`rs-projector-average`), `h_fr = ρ^{2s−6}` (for `L_{2−s}` in the preconditioner). Each product is
`2 Xᵀ (diag(a)⁻¹ diag(H a) − H) X v` (§4.3) with `X = U` or `V`; `H a` precomputed per tree with
the same compression. The far field is rank one per block (App. B, 0th-order multipole); the
near field exact. Vector inputs as 3 columns. Exactness at `χ = 0` against the dense operator is
a test.

### 8.4 Back ends and threads

`Backend::{Dense, Hierarchical { theta, chi, leaf_size, gmres }, Auto { dense_max_faces }}`.
Threads: `FlowOptions::threads` (explicit pool); results bitwise identical for any thread count
(a test at 1, 3, 8 threads).

## 9. Public API

### 9.1 Rust (`ddg_repulsive_surfaces`)

```rust
pub struct Surface { /* validated topology: edges, fans, components, boundary loops, orientation */ }
impl Surface {
    pub fn new(n_vertices: usize, faces: Vec<[u32; 3]>) -> Result<Self, RsError>;
    pub fn validate_positions(&self, positions: &[[f64; 3]]) -> Result<(), RsError>;
    pub fn components(&self) -> &[Component]; pub fn orientation(&self, c: usize) -> &Orientation;
    pub fn boundary_loops(&self) -> &[Vec<u32>]; /* faces(), edges(), n_* ... */
}
pub fn read_obj(path) -> Result<(Surface, Vec<[f64; 3]>), RsError>;   // triangles only (ddg_core::io)
pub fn write_obj(path, &Surface, &[[f64; 3]]) -> Result<(), RsError>;

pub struct Exponents { /* α, β */ }                     // Exponents::new validates (§2.2); Default (6, 12)
pub enum Energy { TangentPoint, Coulomb { alpha: f64 }, Willmore }
pub struct BhOptions { pub theta: f64, pub leaf_size: usize }
pub struct GmresOptions { pub restart: usize, pub tolerance: f64, pub max_iterations: usize, pub side: GmresSide }
pub enum Backend { Dense, Hierarchical { bh: BhOptions, chi: f64, gmres: GmresOptions, eps_rel: f64 },
                   Auto { dense_max_faces: usize } }
pub enum Metric { Sobolev { low_order: bool }, L2(L2Kind), H1, H2 }

pub enum Constraint {
    Barycenter { component: Option<u32>, target: Option<[f64; 3]> },
    Pins { vertices: Vec<u32>, targets: Option<Vec<[f64; 3]>> },
    Area { component: Option<u32>, target: Option<f64> },
    Volume { component: Option<u32>, target: Option<f64> },
}
pub enum BarycenterMode { PerComponent, Global }
pub struct Schedule { pub constraint: usize, pub factor: f64, pub offset: f64, pub steps: usize }
pub struct PinSchedule { pub constraint: usize, pub displacement: Vec<[f64; 3]>, pub steps: usize }
pub enum Penalty { /* §5.2 */ }
pub enum ObstacleKind { TangentPoint, Coulomb { power: f64 } }
pub enum Quadrature { DualArea, PerVertex }
pub struct Cylinder { /* ImplicitSurface */ }

pub enum ArmijoRule { Standard { c: f64 }, Reference { c: f64 } }
pub enum StartRule { Reference, MaxDisplacement { lengths: f64 }, GrowPrevious { factor: f64, first: Box<StartRule> } }
pub struct LineSearch { pub armijo: ArmijoRule, pub start: StartRule, pub tau_min_rel: f64, pub ccd: Option<CcdOptions> }
pub enum NewtonSteps { Fixed(u32), UntilTolerance { max: u32 } }
pub struct RemeshOptions { /* §10 */ }
pub enum Method { ProjectedGradient, Aqp { kappa: f64 }, H1Lbfgs { memory: usize, variant: LbfgsVariant },
                  Bcqn { memory: usize } }                      // §14.2; the last three need a soft problem
pub struct FlowOptions {
    pub exponents: Exponents, pub energy: Energy, pub tp_weight: f64,
    pub metric: Metric, pub backend: Backend, pub method: Method,
    pub constraints: Vec<Constraint>, pub barycenter_mode: BarycenterMode, pub barycenter_shift: bool,
    pub schedules: Vec<Schedule>, pub pin_schedules: Vec<PinSchedule>,
    pub penalties: Vec<Penalty>, pub line_search: LineSearch, pub newton_steps: NewtonSteps,
    pub remesh: Option<RemeshOptions>, pub stop: StopRule, pub max_iterations: usize,
    pub time_limit: Option<Duration>, pub stuck_limit: usize, pub threads: Option<usize>,
    pub reference_compat: bool,        // sets every "reference mode" of §11 at once (comparisons only)
}
// `energy` = ℰ̂ (or Ẽ) without penalties; `penalties[i]` = w_i P_i; `differential` = d(energy + Σ penalties).
pub struct EnergyReport { pub energy: f64, pub penalties: Vec<f64>, pub differential: Vec<[f64; 3]>,
                          pub subcritical: bool, pub pairs_exact: u64, pub nodes_lumped: u64, pub seconds: f64 }
pub fn energy(surface, positions, exponents, eval: &EnergyEval, penalties, threads) -> Result<EnergyReport, RsError>;
pub fn dense_metric(surface, positions, exponents, metric) -> Result<DenseMetric, RsError>;  // B, B0, L_τ, A (tests, oracle)

pub struct Session { /* surface, positions, per-vertex data, targets, history */ }
impl Session {
    pub fn new(surface: Surface, positions: Vec<[f64; 3]>, options: FlowOptions) -> Result<Session, RsError>;
    pub fn step(&mut self) -> Result<StepReport, RsError>;
    pub fn run(&mut self, progress: Option<&mut dyn FnMut(&StepReport) -> Control>) -> Result<FlowReport, RsError>;
    pub fn surface(&self) -> &Surface; pub fn positions(&self) -> &[[f64; 3]];   // change with remeshing
    pub fn pinned(&self) -> &[u32]; pub fn cancel_flag(&self) -> Arc<AtomicBool>;
    pub fn status(&self) -> Status; pub fn options(&self) -> &FlowOptions;
}
pub fn flow(surface, positions, options, progress) -> Result<FlowResult, RsError>;

pub struct StepReport {
    pub iteration: usize, pub status: StepStatus, pub energy_before: f64, pub energy_after: f64,
    pub penalties_after: Vec<f64>, pub gradient_norm: f64, pub sobolev_norm_sq: f64,
    pub tau: f64, pub tau_ccd: Option<f64>, pub halvings: u32,
    pub newton_steps: u32, pub max_violation: f64, pub dropped_rows: Vec<usize>,
    pub linear: LinearReport,          // GMRES iterations per solve, residuals, dense/iterative
    pub remesh: Option<RemeshReport>,  // splits, collapses, flips, refused moves, |V|, |F|
    pub faces: usize, pub seconds: Timings,  // bvh, energy, gradient, factor, solve, line search, ccd, projection, remesh
}
```

**Errors** (`RsError`, hand-written `Display`, `#[non_exhaustive]`): the mesh errors of §2.1,
`PositionCountMismatch`, `NonFinitePosition`, `DegenerateFace { face }`, `NonFiniteEnergy`,
`InvalidExponents`, `InvalidOption { name, value, requirement }`, `TranslationNotFixed
{ component }`, `VolumeNeedsOrientation { component, edge }`, `InfeasibleConstraint`,
`ConstraintOutOfRange`, `MethodNeedsSoftProblem`, `SingularSaddle`, `SolverDiverged { residual }`,
`NonFiniteEnergy { faces: [u32; 2] }` (the first non-finite pair, or the largest term if the sum
overflows), `PolygonFace { face, corners }` (OBJ), `Parse { path, line, message }`,
`Unsupported { what }` (a feature of a later work package, during development only),
`ObstacleError`, `Io`, `Parse { line }`, `Curves(RcError)` (errors from the curves crate's generic
parts, wrapped), `Cancelled`.

### 9.2 Python (`crates/ddg-py`, module `ddg_tools`)

`RepulsiveSurfaceSession(positions, faces, **options)` with `step()`, `run(steps)`, `snapshot()`
(positions and faces, which change with remeshing), `cancel()`; `repulsive_surface_flow(...)`,
`surface_energy(...)`, `repulsive_surface_defaults()`. Options as keywords, mirroring §9.1 and the
curves binding's conventions (buffers in, `array.array` + dicts out, GIL released in
`step`/`run`, unknown option is an error). Typing stubs, README section, tests.

### 9.3 Blender (user decision 2026-11-01)

As for curves: `object.ddg_repulsive_surface_flow` (N steps, one undo step) and
`object.ddg_repulsive_surface_live` (modal, worker thread, Accept/Cancel), panel in the DDG
sidebar tab. Input: the active mesh object's evaluated triangles (loop triangles). Pins from a
vertex group; remeshing never moves or removes pinned vertices (§10). Obstacles: mesh objects of a
collection (tangent-point, weight). Implicit obstacles and attractors from empties (sphere,
cylinder, torus, plane by the empty's display type and scale; repel/attract, power, weight).
Constraints: barycenter (with shift), area and volume (with growth factor and steps), penalties of
§5.2 with weights. When remeshing is on or the topology changes, the mesh data is **replaced** by
the result's triangles; UVs, colours and other attributes are dropped with a warning in the report
(no transfer, user decision); material index and smooth flag are taken from the parent face when
the remesher reports one. With remeshing off and no topology change, positions are written into
the existing mesh and polygons are kept (as `ddg-developability` D24).

## 10. Remeshing (Sec. 6; D16, D17)

A mutable face-slot mesh with per-vertex face lists (a half-edge is a `(face, side)` pair; `remesh` module), built from `Surface` + positions, returning a new
`Surface`, positions, a vertex map (new vertex → old vertices and weights; for per-vertex data)
and a face parent (new face → old face, for Blender's material index). `L₀` = mean edge length of
the **initial** mesh, fixed for the run (paper, reference).

1. **Edge lengths**: split every edge longer than `3L₀/2` at its midpoint (edges existing at the
   start of the pass), then collapse edges shorter than `L₀/2` to their midpoint (a boundary
   endpoint keeps its position; an edge with a pinned endpoint collapses to the pinned vertex;
   both pinned: no collapse), unless the collapse fails the link condition, would create a
   non-manifold vertex or a triangle with an empty interior, flips a face normal by more than
   `π − 2` rad (the reference's foldover threshold, `ref: remeshing.cpp:144–187`), or brings a new
   triangle closer than `gap` to a non-adjacent triangle (only when CCD is on). The paper checks
   foldovers for both splits and collapses; a midpoint split cannot fold a face, so only collapses
   are checked (as the reference).
2. `N = 5` times: **Delaunay flips** (interior edges whose opposite angles sum to more than `π`,
   in `f64`, refused if the flip would make a non-manifold edge or turn a normal by more than
   `π − 2`; a queue with a cap of `100 |V|` flips, reported if hit), then **tangential smoothing**
   `u_i = ρ Σ_{S∋i} a_S (c_S − f_i) / Σ_{S∋i} a_S` with `ρ = 0.5`, `c_S` the circumcenter (the
   barycenter for faces with a boundary edge, as the reference), projected onto the tangent plane
   of the angle-weighted vertex normal, Jacobi update. Pinned and boundary vertices do not move. A
   move that flips a face normal by more than `π − 2` is halved up to 3 times, then refused
   (reported). The paper's order (flip, then smooth); the reference smooths first (reference mode
   available).
**Single-operation conventions** (fixed by the oracle's `remesh-ops`, `Remesh.lean`; the Rust
operations follow them so the fixtures replay exactly): a split appends the midpoint as vertex
`|V|` and each cut face's new face is appended in increasing parent order; a collapse of `(a, b)`
keeps `min(a, b)`, removes `max(a, b)` and shifts higher indices down; its position is the pinned
endpoint, else the single boundary endpoint, else the midpoint (both boundary: midpoint; both
pinned: refused); a flip of edge `p q` with faces `F1 < F2`, `F1` traversing `p → q` and the
opposite vertices `c` (in `F1`) and `d` (in `F2`), gives `F1 := (p, d, c)`, `F2 := (d, q, c)`.
Always checked: the link condition (an interior edge between two boundary vertices fails), no
repeated face, a valid manifold afterwards; a flip that would create an existing edge `c d` is
refused. On a non-orientable surface the two faces of a flip may traverse the edge in the same
direction: the Rust remesher reorders one of them locally first (the oracle refuses such a flip;
fixtures avoid it). The full pass of steps 1–2 may renumber freely (only the vertex map and face
parents are contractual).

3. Per-vertex data through the vertex map: pins stay pinned (never removed; a split edge between
   two pins gives an unpinned midpoint); squared-error targets and attractor membership
   interpolate (the reference does not update the targets, D17).

**As built (WP4).** The foldover test measures how far the **same** face's normal turns under the
operation (`> π − 2` refuses), as the oracle does; the reference instead measures the dihedral angle
to the face across each link edge. A `DegenerateFace` refusal (a zero-area result) is always on for
split, collapse and flip (the stress test produced zero-area faces without it); the report counts
refused splits too. Validation requires `split > 2 · collapse`. A failed revalidation of the
output is `RsError::Internal { what, detail }` (never observed; WP5 adds the variant). Measured: a
stretched icosphere(3) scaled (4, 1, 0.5) goes from 43.5 % of edges outside `[L₀/2, 3L₀/2]` and a
9.2° minimum angle to 0 % and 25.3° after one pass (32.3° after 8); one pass on 40–55k faces takes
0.13–0.20 s on one thread.

`RemeshOptions { split: 1.5, collapse: 0.5, rounds: 5, rho: 0.5, flip_cap_per_vertex: 100,
order: FlipThenSmooth, enabled_ops }`; `remesh: None` disables it (Fig. 7 "without remeshing", and
Blender's "keep topology"). Any remeshing output is revalidated as a `Surface`; a failure is an
internal error with the step's data (never silently repaired).

## 11. Decisions: the paper as printed, the reference, and what we do

| D | topic | paper | reference | ours |
|---|---|---|---|---|
| D1 | mesh type | — | geometry-central (non-orientable OK) | own `Surface`, orientation optional |
| D2 | `D_f` edges | "unit edge vectors" | plain edges | plain edges (misprint) |
| D3 | exponents | `p`, `β = 2p` | free `(α, β)` | free, validated `s ∈ (1,2)`, default (6, 12) |
| D4 | BH gradient | forward terms only (4.2) | exact gradient of `Ẽ` | exact gradient of `Ẽ` |
| D5 | BH admissibility | per face `S` | per leaf, AABB distance | per leaf (reference) |
| D6 | metric weight | — | × TP weight | × `tp_weight` |
| D7 | Willmore H² term | in `A` | unweighted, GMRES only | weighted, in `A` only (in the preconditioner it diverged, WP6) |
| D8 | barycenter | per component + mean-motion shift (Fig. 8), "copies of I" | global, weights `a_v/A`, ε = 1e-2 | `Auto`: global with `B⁰` on, per component + shift with `low_order: false`; weights `a_v/A_c` (§5.1) |
| D9 | volume gradient | "∝ Σ aN" | `Σ aN` (3×) | exact (`⅓ Σ aN` on closed fans) |
| D10 | redundant rows | — | — | rank check, drop + report |
| D11 | growth schedule | "gradually" | after projection | before the step's projection; reference option |
| D12 | Newton correction | 1, repeat if needed | exactly 1 | until tolerance, ≤ 5; `Fixed(1)` option |
| D13 | GMRES | GMRES, x₀ = 0 | Eigen, left, tol 1e-4, restart 30, max 2n | right preconditioned, same constants; left option |
| D14 | preconditioner `ε` | none | 1e-8 / 1e-2 absolute | `1e-10 · tr W / tr M` |
| D15 | Armijo | "standard" | no `‖g‖`, `τ₀ = 2/‖x‖` | standard; start by measurement; reference option |
| D16 | remesh order/checks | flip → smooth; foldover on split/collapse | smooth → flip; collapse only; float angles; flip cap | flip → smooth; collapse + flip + smoothing checks; f64; cap reported |
| D17 | per-vertex data | pins: remeshing off | not transferred; pins break | transferred; pins never touched |
| D18 | exact mode | θ = 0 | broken in the flow (stale BVH, no metric) | `Backend::Dense` works; the reference's θ = 0 is not used as ground truth |
| D19 | `disable_near_field` | drop low-order term | drops all of `B⁰` | `low_order: false` (same) |
| D20 | CCD | suggested (Sec. 10) | none | implemented; default by measurement |
| D21 | implicit quadrature | per point | per vertex, unweighted | `DualArea` default; `PerVertex` option |
| D22 | attractor set | — | UV `u > 1e-10` | explicit vertex set |
| D23 | mesh obstacle | Coulomb-like | two-sided TP + metric diagonal (last obstacle only) | both; TP default; all obstacles |
| D24 | penalty gradients | — | several inconsistent (§5.2) | exact gradients of the values |
| D25 | L² metric | mass matrix | identity | mass; identity option |
| D26 | H¹ L-BFGS | `H₀ = Δ⁻¹`, inner products in Δ | non-standard (dual vectors treated as primal) | standard two-loop with `H₀ = W′⁻¹`; reference variant option |
| D27 | BCQN | Zhu et al. 2018 without barriers | Frobenius norm, `B = A^{4/3}` | Zhu's blend with a spectral-norm estimate (power iteration), no barrier; reference variant option |
| D28 | AQP | Kovalsky et al. 2016 | κ = 100 fixed, current Laplacian, no safeguards | κ = 100, current Laplacian (remeshing changes the rest pose), momentum reset on energy increase and remeshing |
| D29 | stop rule | — | iteration/time limit only | + energy stall (§7.4) |
| D30 | shared code | — | — | depend on `ddg-repulsive-curves` (its SPEC §16 had proposed extracting a shared crate; the generic parts are public and stable, so we use them in place and do not edit that verified crate) |
| D33 | preconditioner form | `L′⁻¹ L_{2−s} L′⁻¹` | same, multipliers zeroed | augmented `L′⁻¹ [[L_{2−s}, C_Qᵀ],[C_Q, 0]] L′⁻¹` (§6.3); printed form as reference mode |
| D32 | preconditioner's `L_{2−s}` | face form | face form, left GMRES | vertex pairs (§6.3); face form as reference mode |
| D31 | Fig. 12 meshes | not published | not in the repo | generated analogues at the paper's face counts (§14.4), recorded as such |

`reference_compat: true` sets D8 `Global` + ε 1e-2, D9 3×, D11, D12 `Fixed(1)`, D13 left, D14
absolute, D15 reference rule and start, D16 smooth→flip, D21 `PerVertex`, D25 identity, D26/D27
reference variants, and the penalty value conventions of the scene keywords — used only for the
black-box comparison (§14.7).

## 12. Formula ids

Rust functions carry `[formula: <id>]`; each id is a Lean theorem in
`lean/proofs/DdgProofs/RepulsiveSurfaces/*.lean` ("proved") or a heading in `lean/proofs/CITED.md`
("cited").

| id | statement | Rust (planned) | Lean / CITED |
|---|---|---|---|
| `rs-exponents` | `s = (β−2)/α`; `= 2 − 2/p` for `β = 2p`; `2σ+2 = 2s`; `2(2−s)+2 = 6−2s`; `s ∈ (1,2) ⇔ α+2 < β < 2α+2` | `Exponents` | proved, `Basic.lean` |
| `rs-exponent-range` | finite energy ⇔ `W^{s,p}` embedded (Strzelecki–von der Mosel 2013, Blatt 2013); barrier for `p > 2n` | `Exponents` (doc) | cited |
| `rs-differential-order` | `dℰ` has order `2s`; `B` matches it (Sec. 2.3, 2.4) | docs | cited |
| `rs-face-quantities` | `a = ½|n|`, `X`, `N`; `a_v` sums to the total area | `geometry` | proved, `Basic.lean` |
| `rs-face-area-gradient` | `∂a/∂f_i = ½ N × e_jk` | `geometry` | proved, `Geometry.lean` |
| `rs-face-normal-variation` | `∂N/∂f_i = (1/2a)(e_jk × N) Nᵀ` | `geometry` | proved, `Geometry.lean` |
| `rs-kernel` | Eq. 16; even in `N`; rigid invariance; homogeneity `c^{α−β}`; `|Pv|² = (N·v)²` for `P = NNᵀ` | `kernel` | proved, `Kernel.lean` |
| `rs-kernel-gradient` | §3.2 kernel partials | `kernel` | proved, `Kernel.lean` (if rpow calculus is heavy: FD + oracle, no tag) |
| `rs-energy` | Eq. 10 as an ordered double sum, `S ≠ T` | `energy::exact` | proved (definition), `Energy.lean` |
| `rs-energy-gradient` | chain rule of §3.2 | `energy::exact` | proved, `Energy.lean` |
| `rs-energy-scaling` | `ℰ̂(cf) = c^{4+α−β} ℰ̂(f)` | `energy` (doc) | proved, `Energy.lean` |
| `rs-gradient-identities` | `Σ g = 0`, `Σ f × g = 0`, Euler identity | `energy` (doc) | proved, `Energy.lean` |
| `rs-bh-energy` | Eq. 17 with the leaf rule; `θ = 0` ⇒ exact | `energy::barnes_hut` | proved, `BarnesHut.lean` |
| `rs-bh-aggregate-adjoint` | §3.4 percolation = gradient of `Ẽ` | `energy::barnes_hut` | proved, `BarnesHut.lean` |
| `rs-box-distance` | `d² = Σ max(0, lo−hi)²` is a lower bound of the set distance | `spatial` | proved, `BarnesHut.lean` |
| `rs-df-gradient` | §4.1 `D_f` = P1 gradient; unit-edge counterexample | `metric` | proved, `Metric.lean` |
| `rs-sobolev-forms` | Eq. 11–13 increments = the quadratic forms; PSD | `metric::dense` | proved, `Metric.lean` |
| `rs-low-order-symmetrisation` | ordered sum = symmetric kernel | `metric` | proved, `Metric.lean` |
| `rs-kernel-matrix-form` | App. A identity (§4.3), `U` and `V` | `metric`, `hierarchical` | proved, `Metric.lean` |
| `rs-hierarchical-null-space` | `L·1 = 0` with any symmetric compressed `H` when `diag` uses the same `H` | `hierarchical` | proved, `Metric.lean` |
| `rs-projector-average` | `vᵀP̄v = Σ w_S |P_S v|²` for a convex combination | `hierarchical` | proved, `Metric.lean` |
| `rs-fractional-operator` | Eq. 5/6 fractional Laplacian | docs | cited |
| `rs-preconditioner-order` | `W⁻¹(M Δ^{2−s})W⁻¹ = Δ^{−s} M⁻¹` for `Δ = M⁻¹W` (commuting powers of a diagonalizable `Δ`) | `precond` | proved, `Precond.lean` (matrix identity; the approximation of `B` cited) |
| `rs-fractional-preconditioner` | augmented `P⁻¹ = L′⁻¹ N L′⁻¹`, `N = [[L_{2−s}, C_Qᵀ],[C_Q, 0]]` (D33): equals `L′⁻¹ + L′⁻¹ blockdiag(L_{2−s} − K, 0) L′⁻¹`, symmetric, nonsingular iff `N` is, `C_Q (P⁻¹[r; b])_x = b` | `solve::precond` | proved, `Precond.lean` (`RS.aug_precond_props`; the printed form's `frac_precond_props` kept untagged) |
| `rs-preconditioner-approx` | `Ã⁻¹` approximates the inverse of the high-order term (Sec. 5.2) | docs | cited |
| `rs-schur-direction` | Eq. 19–20: `x` solves the saddle, `Cx = 0` | `solve::schur` | proved, `Saddle.lean` |
| `rs-schur-correction` | Eq. 21: `C h = Φ` | `solve::schur` | proved, `Saddle.lean` |
| `rs-nested-saddle` | positional rows inside `A′` + Schur rows = the full saddle (Sec. 7.2) | `solve::schur` | proved, `Saddle.lean` |
| `rs-barycenter` | `Φ_C`, its frozen-weight Jacobian, translation projection exact | `constraints` | proved, `Constraints.lean` |
| `rs-area-gradient` | `∇A` | `constraints` | proved, `Constraints.lean` |
| `rs-volume-gradient` | `∇V = ⅙ Σ f_j × f_k`; `= ⅓ Σ a_S N_S` on a closed fan (the reference's 3×) | `constraints` | proved, `Constraints.lean` |
| `rs-penalties` | values and gradients of §5.2 (area/volume deviation, powers, squared error, boundary length, turning angle, implicit potentials) | `penalties` | proved, `Penalties.lean` |
| `rs-willmore` | `fᵀ W M⁻¹ W f` and its gradient | `penalties::willmore` | proved (structure) / FD, `Penalties.lean` |
| `rs-armijo-reference` | the reference's test = standard Armijo with `c/‖g‖` | `step::line_search` | proved, `Step.lean` |
| `rs-circumcenter` | barycentric weights `a²(b²+c²−a²) : …` give the circumcenter | `remesh` | proved, `Remesh.lean` |
| `rs-smoothing-weights` | `Σ_S a_S = 3 a_v`: the reference's `½ Σ a_S(c−x)/(3a_v)` is `ρ = ½` | `remesh` | proved, `Remesh.lean` |
| `rs-delaunay-test` | opposite-angle sum `> π` ⇔ cotan weight `< 0` | `remesh` | proved, `Remesh.lean` |
| `rs-ccd-advancement` | conservative advancement with moving triangles is a lower bound of the first contact time | `ccd` | proved, `Ccd.lean` (extends `rc-ccd-advancement`) |
| `rs-lbfgs-two-loop` | two-loop recursion with `H₀ = W′⁻¹` applies the BFGS inverse | `optimize` | proved, `Optimize.lean` (or cited: Nocedal–Wright) |
| `rs-aqp` / `rs-bcqn` | Kovalsky 2016 momentum; Zhu 2018 blend | `optimize` | cited |
| `rs-energy-consistency` | `ℰ̂(f_h) → ℰ(f)` empirically between O(h) and O(h²) (Fig. 10) | docs | cited |
| `rs-chen-holst` | the remeshing scheme of Chen–Holst 2011 | docs | cited |

## 13. Verification

### 13.1 Layers (`CLAUDE.md`, adapted)

No unique minimizer to certify (as for the cut flow and curves): (1) analytic values (the energy
of a round sphere has a closed form for the smooth surface; the discrete energy of two parallel
squares; the Df of linear functions), (2) finite differences of every energy, penalty and
constraint, (3) invariances and scaling (§3.3), (4) dense vs hierarchical operators and energies
at `θ = χ = 0` (exact equality to round-off) and their error decay in `θ` (Fig. 10 right:
`O(θ²)`), (5) the Lean oracle, (6) the Lean proofs, (7) the reference's own program as a black box
(§14.7).

### 13.2 Finite differences

Central differences with step `1e-6 ℓ̄`, relative error `< 1e-6` on random perturbed meshes
(icosphere, torus, a mesh with boundary, a non-orientable band, two components), for: the exact
energy, the BH energy (its gradient against its own FD with frozen lists), every penalty and
constraint value, the Willmore term.

### 13.3 Metric and solves

`A` symmetric PSD, kernel = per-component translations (and nothing else on connected meshes,
smallest nonzero eigenvalue reported); hierarchical vs dense product at `χ = 0`: `1e-12`
relative; at `χ = 0.5`: `< 5 %` (reported); preconditioned GMRES reaches `1e-4` in fewer than 50
iterations on the Fig. 13 set (measured, §14.5); dense Schur = full saddle LU to `1e-10`; Newton
correction reduces `|Φ|` quadratically on area/volume.

### 13.4 Flow

Energy decrease every accepted step; constraints held to tolerance; a torus flows to a torus of
revolution (Fig. 15 g = 1), a sphere stays round, a perturbed sphere becomes round; handcuffs do
not intersect (checked by an exact triangle–triangle intersection test on the result: this is a
test utility, not a library guarantee unless CCD is on); bitwise equality across thread counts;
session step-by-step = one-call flow; remeshing keeps manifoldness and pins.

### 13.5 Lean oracle (`lean/oracle`, no Mathlib; `DdgOracle/RepulsiveSurfaces*.lean`)

Independent Float implementation from this file only (shares nothing with the Rust code). Tasks:
`energy` (exact, value + gradient), `metric` (dense `L_τ`, `B`, `B⁰`, `A`), `df`, `constraints`
(values + Jacobians), `saddle` (dense nested saddle + Schur direction + correction), `precond`
(dense `L′` and `P⁻¹ r`), `step` (one dense step with standard Armijo, no remeshing), `remesh-ops`
(one split/collapse/flip on a small mesh, positions + connectivity). Golden fixtures in
`lean/oracle/fixtures/rs/`, `MANIFEST.tsv`; Rust test `tests/oracle.rs` replays them always and
fuzzes against the executable when asked (`#[ignore]`), as `rc`.

### 13.6 Lean proofs and citations

`lean/proofs/DdgProofs/RepulsiveSurfaces/{Basic,Geometry,Kernel,Energy,BarnesHut,Metric,Precond,
Saddle,Constraints,Penalties,Step,Remesh,Ccd,Optimize}.lean`; `CITED.md` section
`## Repulsive Surfaces`; `scripts/check-formulas.sh` must pass; main theorems added to the axiom
audit. Lean workspace rules: `cd lean/proofs` (never `lake --dir`), never `lake update`, never
rebuild Mathlib (worktrees symlink `.lake/packages`).

## 14. Measurement plan

### 14.1 Consistency (Fig. 10)

The perturbed torus of Sec. 8.1 (`f₀ + Φ∘f₀`, a fixed seeded trigonometric polynomial). The
"true" energy by our own adaptive cubature over `[0, 2π]⁴` (the integrand is bounded on a smooth
embedded surface; tensor Gauss–Legendre with subdivision of the diagonal cells, parallel; two
independent rules must agree to 6 digits) — no Mathematica. Discrete energies at 6+ resolutions
(sampled on a squeezed Delaunay triangulation of the parameter domain, remeshed and projected
back), and the Barnes–Hut error against `θ`. Expected: rate between `O(h)` and `O(h²)`, and
`O(θ²)`.

**Results of record (WP16a, `examples/consistency.rs`, `NOTES-WP16A.md`).** Perturbed torus
`R = 1`, `r = 1/3`, `PerturbedTorus(degree 2, amplitude 0.2, seed 16)` (embedded, area 13.1631).
Smooth energy by our own cubature over `[0, 2π]⁴` (offset coordinates, trapezoid in the base point,
Duffy and geometric-cell rules for the diagonal, which agree to 4e-8): `ℰ(6,12) = 460.71622`,
`ℰ(4,8) = 96.400663`; torus of revolution `ℰ(6,12) = 103.770432`, `ℰ(4,8) = 61.645093`. Discrete
energy on squeezed grids remeshed 4 passes with projection, 648–85,554 faces (`h` 0.29–0.033):
relative error 2.0 → 2.7e-2, fitted rate 1.97 (1.86 over the last four points), torus 2.05, (4, 8)
1.97 — at the `O(h²)` end of the paper's "between `O(h)` and `O(h²)`". Barnes–Hut (leaf 8): error
`~ θ²` for `θ ≥ 0.35` (fits 1.9–2.4; at `θ = 0.5`: 1e-4 / 4.8e-3 / 8.2e-3 on 11k / 43k / 86k
faces); for `θ ≤ 0.25` the signed error changes sign and is not monotone (`|error| ≤ 2e-4` at
`θ = 0.125`) — cancellation between over- and under-estimated far fields, not investigated
further.

### 14.2 Comparisons (Fig. 5, 12; Sec. 8.2)

Methods: `Hs` (ours), L², H¹, H² projected gradient, AQP, H¹ L-BFGS, BCQN; all with the same
energy/gradient code and remeshing; hard barycenter + area constraints for the projected methods,
the area penalty (`AreaDeviation`, weight 100 as the reference's `scene_soft.txt`) for AQP and the
L-BFGS methods (Sec. 8.2.1). Fig. 5: one surface at ~11k / 40k / 135k faces, 300 steps each, all
methods. Fig. 12: the five meshes (§14.4 analogues), 2400 s (3600 s for the two tunnels) per run,
one thread, energy against wall time, exact energy for the plots. Runs in parallel (one per core,
memory permitting) with one log per run; a timed subset alone afterwards (as the curves
evaluation, user decisions Q4/Q7 there).

**Results of record (WP16b, `examples/compare_methods.rs`, `scripts/rs-eval.sh`,
`NOTES-WP16B.md`).** One thread per run, 56 runs (35 Fig. 12 + 21 Fig. 5) in parallel; the
machine's load average was 79 before the runs started and 51 at the end (outside jobs), so absolute
times are inflated and per-second comparisons favour the methods with cheaper steps. Energies are
reported as `E · A/A₀` (exact energy at the end), because the soft-penalty methods let the area grow
up to 2.2×. Fig. 12 end values, `Hs` vs the best other method: linked handcuffs 7.7 vs 26.0,
doubly-wound 7.2 vs 60.0, 6×-wound 902 vs 1179 (BCQN), figure-8 tunnel 353 vs 341 (BCQN), trefoil
tunnel 158 vs 154 (H¹ L-BFGS). `Hs` is lowest **per iteration** on all five meshes and lowest **per
second** on the first three; on the two tunnels H¹ L-BFGS / BCQN are as good or slightly better per
second (`Hs` 20–36 s/step on 24–54k faces under that load, 20–30 GMRES iterations per solve on meshes
3–5). Fig. 5: `Hs` progress at step 100 is resolution-independent (41.1 / 39.5 / 39.1 at 11k / 40k /
135k faces); L² slows 2.6×; L-BFGS/BCQN slow mildly; H¹, H² and AQP are flat on this easy genus-2
surface, so that part of the paper's claim is not reproduced. No run ended Stuck or with an error;
AQP's energy spiked to ≈ 6e6 on mesh 3 (reproducer in the notes; the paper also reports AQP
blow-ups). The paper's wall-clock claim is therefore only partly reproduced under these conditions;
a rerun of the two tunnels with fewer concurrent jobs would separate method from machine load.

### 14.3 Defaults confirmed by measurement

Line-search start rule (§7.2), CCD default (§7.3), `dense_max_faces`, GMRES side and tolerance,
`ε_rel`, Newton cap, Willmore term in the preconditioner, `θ`/`χ` — each with a table in the
README (criterion: robustness, correctness, speed, in that order).

### 14.4 Applications and test meshes (Sec. 9; user: all)

Generators (`generate` module + `examples/`): genus-g surfaces (linear handle arrangement),
tubes around curve networks (implicit union of capsules/smooth tubes, meshed by marching
tetrahedra, then remeshed) for: handcuffs (linked, doubly- and 6×-wound), figure-8 and trefoil
tunnels, the knotted genus-2 table of Ishii et al. 2012 (Fig. 16; at least 4₁, 5₁–5₄, 6₁–6₁₆ as far
as the table's diagrams can be transcribed; graphs optionally relaxed by
`ddg-repulsive-curves` first), the "H" text of Fig. 22, punctured torus (Fig. 20), cobblestone
spheres in a slab (Fig. 25), a polygon soup with holes (Fig. 23, from a mesh with faces deleted).
Reference scenes (handcuffs-and-pole, double torus, walnut, Klein bottle, bunny/spot point clouds,
borromean, Hopf, …) are imported by a scene-file reader (`scene` module, the reference's
`FORMAT.txt` grammar incl. the undocumented keywords) after a fetch script copies them into
`data/repulsive-surfaces/` (gitignored; `DATASETS.md` rows). Applications: Fig. 15 (g = 0…11),
16, 17 (attracting plane), 18 (handle off a pole; implicit cylinder, free barycenter), 20
(eversion: zero volume + fixed area + boundary penalties), 21–22 (pins, squared error, volume
growth), 23 (shrink-wrap: point cloud + soup), 24 (nested envelopes: volume 1.15× per level with
the previous level as obstacle), 25 (walnut, cobblestones), 26 (p < 4, `low_order: false`; Klein
bottle failure reproduced). Each: an example, a smoke test at reduced size, and a showcase entry
(`showcase/repulsive-surfaces/`).

### 14.5 Performance (Fig. 13; the floor)

Per-step time at 4 threads on the 14 Fig. 13 models (or analogues at the same face counts): at
least as fast as the paper (≈ linear, ≈ 2.5 s/step at 45k faces on a Threadripper 3990X); plus
thread scaling (1–64) and a breakdown (BVH, energy, gradient, factor, GMRES products, line search,
remesh). Benches are plain `harness = false` mains printing CSV, run alone.

### 14.6 Blender end-to-end

Headless tests of both operators (N steps with and without remeshing, pins kept, obstacles,
implicit empties, undo, cancel restores), and a windowed live test on Xvfb, as for curves.

### 14.7 Black-box comparison with the authors' program

Build outside the repository (`/var/data/discrete-geometry/RepulsiveSurfaces-build/`): fetch the
three submodules at the pinned commits, MKL from a user-space prefix (micromamba or the `mkl-devel`
wheels via `uv`), a headless driver (~150 lines, polyscope removed) that runs a scene for N steps,
logs energy/time and writes OBJs. Run with `--threads 1` (its squared-error race). Compare with
`reference_compat: true`: dense operators on small meshes (Gram matrices, Jacobians, one direction)
as parity fixtures (`tests/reference_fixtures.rs`), and trajectories/final shapes on the shared
scenes side by side (not matched step for step: remeshing amplifies round-off).

**Built (WP18a, 2026-11-02).** `/var/data/discrete-geometry/RepulsiveSurfaces-build/` (README there):
GCC 15.2, Release `-O3 -march=native`, MKL 2026.0.2 from `uv` wheels with `mkl_gnu_thread` +
libgomp (no libiomp5), geometry-central only (polyscope cloned, not built); submodules
geometry-central `00d09a5`, polyscope `feaf66b`, libgmultigrid `8d01391`, all MIT (licences in
`licenses/`). One patch (`gcc-compat.patch`): the reference compiles only with clang; the patch
removes `aligned(...)` hints from `#pragma omp` lines and a duplicated `P_min` in a `shared()`
clause — no numerics. Driver `build/bin/rsurfaces_headless` (`--scene`, `--iterations`,
`--obj-every`, `--out`, `--threads` (1), `--theta` (0.5), `--no-remesh`, `--dump-initial`
[`--dump-metric`], `--coulomb`, `--mult-alg`; writes `steps.csv`, OBJs, the reference's stdout).
Smoke, 20 steps, 1 thread: EntangledDoubleTorus 1036 ms/step, energy 377.0 → 194.5; HopfLink
609 ms/step, 849.3 → 313.7; Bunny 3501 ms/step, 162452 → 100.06. **The reference cannot load the
Klein bottle meshes of its own `scenes/KleinBottle/`** (its geometry-central version aborts with
"duplicate edge in list" on non-orientable input): Fig. 26 (bottom) is not reproducible with the
published `master`; our `Surface` accepts them (D1). The authors' unmerged branch `nonorientable`
(9 commits ahead of `master`, 45 behind, same geometry-central `00d09a5`) has commit `ce89c2c`
("changes necessary to run nonorientable meshes"): geometry-central's general `SurfaceMesh` instead
of the orientable `ManifoldSurfaceMesh`, and edge splits and collapses switched off (the general mesh
type of that version does not support them; flips and smoothing stay). Fig. 26 (bottom) was most
likely made with it. **Comparison variant (user, 2026-11-02):** a second binary,
`rsurfaces_headless_nonorientable`, built from `master` + `ce89c2c` only (not the whole branch, so the
Klein-bottle change is the only difference from the main binary). On non-orientable meshes the
comparison runs our flow with splits and collapses off too (`RemeshOptions::enabled_ops`), so both
sides remesh alike; a second run of ours with full remeshing shows what the published code cannot
do. Built (2026-11-02): `ce89c2c` applied to the build's sources with one
hand-resolved hunk (the two reader calls in `rsurface_types.h`), the driver compiles unchanged
against `SurfaceMesh`. Smoke, 20 steps, 1 thread, |V| and |F| constant throughout:
`klein_bottle_blocky.obj` (1,800 faces) 559 ms/step, energy 793.91 → 206.31;
`klein_bottle_2.obj` (16,384) 10,991 ms/step, 387.20 → 171.24; `klein_bottle.obj` (16,384): the
energy **rose** at step 1 (309.77 → 323.29) and step 2 ran over 8 minutes at full CPU (killed, not
investigated); EntangledDoubleTorus 377.047 → 194.4899 (main binary 194.4891).

## 15. Work packages

Each package is one agent, on its own worktree branch merged by the orchestrator. "Files" are the
only paths it may edit; WP1 creates stubs of every module so parallel packages never touch the
same file. `src/lib.rs` is edited only by WP1 and the orchestrator. Briefs repeat the standing
rules (CLAUDE.md, parallel independent runs, own `CARGO_TARGET_DIR`, Lean workspace rules, commit
often, notes file).

| WP | content | files | depends on | model |
|---|---|---|---|---|
| WP1 | crate skeleton, `Surface` + validation (orientation classes), geometry (areas, normals, barycenters, `W`, `M`), exponents, kernel, exact energy + gradient, OBJ I/O, errors, module stubs, LICENSES, README stub | `Cargo.toml`, `src/{lib,error,surface,geometry,io,kernel}.rs`, `src/energy/{mod,exact}.rs`, stubs, `tests/{surface,energy}.rs` | — | Opus |
| WP2 | dense metric (`D_f`, `L_τ`, `B`, `B⁰`, `A`, L²/H¹/H²), constraints (Jacobians, rank, per-component), dense nested saddle + Schur + correction, GMRES (generic, right/left), dense preconditioner | `src/metric/{mod,dense}.rs`, `src/constraints/*`, `src/solve/{mod,dense,schur,gmres}.rs`, `tests/{metric,constraints,gmres}.rs` | WP1 | Opus |
| WP3 | Barnes–Hut energy + exact adjoint gradient on `Bvh<3>`, refit, determinism | `src/spatial/*`, `src/energy/barnes_hut.rs`, `tests/barnes_hut.rs` | WP1 | Opus |
| WP4 | remeshing (mutable half-edge, split/collapse/flip/smooth, legality, pins, boundary, vertex map, face parent) | `src/remesh/*`, `tests/remesh.rs` | WP1 | Opus |
| WP5 | session, step, line searches, projection, schedules, shift, stop rules, comparison metrics L²/H¹/H² | `src/step/*`, `tests/flow.rs` | WP2, WP3, WP4 | Opus |
| WP6 | block cluster tree kernels, hierarchical `A` product, sparse preconditioner (`L′` LDLT), iterative Schur, obstacle diagonal, `Backend::Auto` | `src/metric/hierarchical.rs`, `src/solve/{iterative,precond}.rs`, `tests/hierarchical.rs` | WP2, WP3 | Opus |
| WP7 | penalties + obstacles (TP and Coulomb mesh obstacles with BH, point cloud, implicit + `Cylinder`, Willmore + H² term, Coulomb energy) | `src/penalties/*`, `src/obstacle.rs`, `src/implicit.rs`, `tests/penalties.rs` | WP1, WP3 | Sonnet |
| WP8 | CCD (moving triangles, vertex–face, edge–edge, broad phase, obstacles), line-search cap, cost measurement | `src/ccd.rs`, `tests/ccd.rs`, `benches/ccd.rs` | WP5 | Sonnet |
| WP9 | comparison optimizers: AQP, H¹ L-BFGS (both variants), BCQN (both), soft problems | `src/optimize/*`, `tests/optimize.rs` | WP5 | Sonnet |
| WP10 | Lean proofs + CITED | `lean/proofs/DdgProofs/RepulsiveSurfaces/*.lean`, `…/RepulsiveSurfaces.lean`, one import line, `CITED.md` section, audit entries | SPEC | Opus |
| WP11 | Lean oracle + fixtures; Rust replay/fuzz test | `lean/oracle/DdgOracle/RepulsiveSurfaces*.lean`, dispatch lines in `Main.lean`, `lean/oracle/fixtures/rs/`, oracle README row, `tests/oracle.rs` | SPEC; Rust test after WP5 | Sonnet |
| WP12 | generators (genus-g, tubes via marching tetrahedra, text, soup, cobblestones), scene reader, fetch script, `DATASETS.md` rows | `src/generate/*`, `src/scene.rs`, `scripts/rs-fetch-scenes.sh`, `tests/generate.rs` | WP1 (WP7 for scene penalties) | Sonnet |
| WP13 | Python bindings | `crates/ddg-py/src/repulsive_surfaces.rs`, `src/python/repulsive_surfaces.rs`, registration lines, `Cargo.toml` line, `.pyi`, README section, `python-tests/test_repulsive_surfaces.py`, `scripts/py-test.sh` line | WP5, WP7 | Sonnet |
| WP14 | Blender operators | `blender/ddg_uv_unwrap/{repulsive_surface_operator,repulsive_surface_io}.py`, registration lines, `blender/tests/{test,ui}_repulsive_surface.py`, `scripts/blender-test.sh` lines, README section | WP13 | Sonnet |
| WP15 | applications + showcase (§14.4) | `examples/*`, `tests/applications.rs`, `showcase/repulsive-surfaces/` | WP5–WP7, WP12 | Sonnet |
| WP16 | evaluation (§14.1–14.3) | `examples/{consistency,compare}.rs`, `scripts/rs-eval*.sh`, `scripts/rs_eval_*.py`, results in this file and README | WP6, WP9, WP12 | Sonnet |
| WP17 | benches (§14.5) | `benches/*`, `[[bench]]` entries | WP6 | Sonnet |
| WP18 | reference build + driver (outside the repo), parity fixtures, black-box comparison | `tests/reference_fixtures.rs`, `tests/fixtures/reference/`, `COMPARISON.md` | WP5, WP6 | Sonnet |

Rounds (≤ 3 agents; timed work alone at the end):

1. WP1, WP10, WP11 (Lean side).
2. WP2, WP3, WP4.
3. WP5, WP6, WP7.
4. WP8, WP9, WP12.
5. WP13, WP15, WP18 (+ WP11 Rust side).
6. WP14, WP16 phase A; then alone: WP17, WP16 timed phase, §14.3 sweeps that time things.

## 16. Open items

(Filled in as work packages report.)

* ~~Lean: restate `rs-fractional-preconditioner` for the augmented form~~ — done (WP10b,
  `RS.aug_precond_props`, exactly as §6.3 states it).
* The obstacle part of the preconditioner (reference `fr` diagonal) is not implemented (WP6).
* ~~Naming after WP6~~ — done by WP5: `FractionalForm::{FacePairs, VertexPairs}` (default
  `VertexPairs`) everywhere; `DEFAULT_DENSE_MAX_FACES = 1000`; `RsError::Internal` added.
* **As built (WP5), API differences from §9.1:** `FlowOptions::barycenter_mode:
  Option<BarycenterMode>` (`None` = `Auto`); `Backend` as WP6 built it (`solve/iterative.rs`) plus
  `FlowOptions::energy_eval`; added `StartRule::Absolute`, `StopRule::Off`,
  `LineSearch::max_halvings`, `schedule_timing`, `ccd_obstacles: Vec<Arc<Obstacle>>`,
  `Session::{face_origin, constraints}`, `FlowOptions::effective()`; scene runs via
  `scene_options`/`run_scene` (per-step CSV) and `examples/run_scene.rs`. `w_TP` weights both the
  energy and the metric (the reference's reading; the oracle's `metric` task weights only the
  metric — its fixtures use `w_TP = 1`, so they agree). The ⅔ `τ_max` cap of §7.3 applies only when
  CCD finds a contact before `τ₀` (otherwise every start would shrink). `reference_compat` cannot
  reproduce D9 (3× volume row), the absolute `ε` and `τ_min`, or the D24 value conventions; the
  black-box comparison states this per scene.
* **CCD default (D20), measured by WP5:** the median CCD share of a step is 0.5–2.0 % on the five
  smoke scenes, so CCD is **on** by default, per the user's rule. Two effects to watch:
  (1) on EntangledDoubleTorus the run with CCD reaches `F = 42` after 200 steps, without CCD 8.2
  (no intersection at any step end in either). **Resolved (CCD investigation,
  `NOTES-CCDINV.md`): the run without CCD tunnels.** Steps 1–5 are identical; in step 6 it
  accepts `τ = 3.72` (a vertex moves 14.3 `ℓ̄₀`) and an exact swept test finds 2226 crossings in
  that step (488 vertex–triangle, 1738 edge–edge): the small ring, linked through the large ring,
  is pulled straight through the large ring's tube. No other step, line-search trial, projection
  or remeshing pass crosses (remeshing moves vertices ≤ 0.49 `ℓ̄₀`, sheets stay ≥ 2 `ℓ̄₀` apart).
  The end states are very likely in the same isotopy class (the linked-ring start is still a
  standard genus-2 surface); the run with CCD untangles legitimately by step ≈ 270, and at the
  scene's 700 steps `F` = 7.633 (CCD) vs 7.607 (none), same shape. So the paper's algorithm as
  printed (no CCD) can pass a surface through itself on the paper's own kind of example; CCD keeps
  all 200 steps crossing-free at ≈ 0.9 % of a step (17 ms of 1.9 s). Gap 0.02 / 0.01, cap 0.9 or
  no collapse guard change `F` after 200 steps only within trajectory noise (35.6–41.8), so the
  defaults stay (gap 0.05 `ℓ̄₀`, cap ⅔, guard on). (2) an immersed input (the Klein bottle) has touching pairs, which
  stall CCD; `scene_options` turns CCD off for subcritical exponents. WP5 smoke runs (4 threads,
  without CCD): HopfLink 9.2k faces 563 ms/step, 32 GMRES iterations/step; EntangledDoubleTorus
  10.5k 609 ms, 21; KleinBottle 1.8k 147 ms, 46; HandcuffsPole 10.7k 703 ms, 31; Bunny 28.6k
  2357 ms, 34 (GMRES 1.33 s of it); every step remeshed, none stuck.
* Step time (WP6, 45k faces, 4 threads, run alone; barycenter + area + volume = 3 solves): BH build
  + energy/gradient 0.32 s, face BCT build 0.57 s, factorization 0.04 s, vertex tree 0.10 s, GMRES
  1.15 s (29 products at 25 ms) — 2.2 s for the direction against the paper's 2.5 s for a whole
  step (`χ = 1`: 1.6 s). The BCT build and the products are the targets of WP17.
* **CCD on inputs that already intersect** (Klein bottle immersion; Fig. 26): implemented (band
  rule). At the start of each step the faces that intersect are found (exact triangle–triangle
  test over a BVH, `ccd::intersecting_pairs`); `I` = their vertices widened by one ring; a pair whose
  two elements both have a vertex in `I` is unconstrained (`exempt_pairs` in the report), every
  other pair keeps the gap rule; the collapse guard does the same. An embedded mesh has an empty
  band and the old behaviour (bitwise). CCD now stays on for every exponent. Klein bottle, 200 steps:
  `F` = 25.3 with CCD + band, 155.7 without CCD (631 with CCD and no band). Cost of the start test:
  7–17 ms at 10–50k faces on 8 threads (≈ 0.35 of the CCD itself).
* **Black-box comparison (WP18b, `COMPARISON.md`).** Parity: all-pairs energy 2e-16–7e-16, its
  differential ≤ 1.1e-14, the dense constrained metric ≤ 1.1e-14 of its largest entry, barycenter
  and area rows ≤ 2e-16; the volume row differs by exactly 3× (D9); Barnes–Hut at `θ = 0.5` differs
  by 3e-6–3e-5 (different tree splits). Trajectories (11 scenes): `reference_compat` tracks the
  reference's energy curve to ≈ 1 % and the final shape to ≈ 0.5 % (EDT, Hopf, handcuffs, bunny,
  walnut, blocky Klein bottle, lung); our defaults reach the same plateaus sooner (EDT `F` 8.3 vs
  19.9 at 300 steps); Borromean rings: the reference stalls at `F ≈ 7302` from step 50 and its mesh
  self-intersects in steps 1–9, ours reaches ≈ 4540 with a symmetric shape and no intersection.
  The reference's volume lags its target by −4.5e-6 (ours −4e-8). Klein bottles: master cannot load
  any; the non-orientable variant flows the blocky one (ours within 0.2 % with splits/collapses off),
  stalls on `klein_bottle.obj` (ours runs it), and on `klein_bottle_2.obj` reaches 161.7 vs ours
  117.9 (full remeshing). **Speed (1 thread, alone): ours is slower per step** — 1.1–1.8×
  (`reference_compat`), 1.2–2.3× (defaults); Hopf link 2.4–4.4× because we run 59 GMRES iterations
  per step against the reference's ≈ 11 (its left GMRES stops early with a large true residual,
  D32/D33); faster on walnut (0.7×) and the blocky Klein bottle (0.8×). Target of WP17.
* **Defects from WP18b (`NOTES-DEFECTS.md`):** (1) fixed: `reference_compat` with splits/collapses off on
  the blocky Klein bottle aborted at step 300 with `SolverDiverged { residual: 1.466e3 }`; cause: a flip at
  the end of step 299 made needle faces out of an edge the flow had shrunk to 7e-4, on which the printed
  preconditioner's left GMRES fails; a diverged hierarchical solve now falls back to right/augmented GMRES,
  then a dense LU (§6.3, `LinearReport::fallbacks`). (2) Reconstruction: the reference adds the metric
  diagonal `O_diag` for point clouds too (zero obstacle normals); `reference_compat` now does (first step
  `τ` 3.40e-4 against 3.37e-4, objective within 1.5 % to step 125); the final self-energy (4267 against
  14.9 at step 150, concentrated in creases near the cloud) is still unexplained.
* **Blender (WP14) notes:** `Mesh.clear_geometry()` in Blender 5 also deletes vertex groups (the
  pin group is rebuilt, others restored on Cancel); CCD auto-off on intersecting input used
  Blender's `BVHTree.overlap` before the band rule existed — with the band rule CCD can stay on;
  a free flow without an area constraint inflates the surface (the energy decreases with scale) —
  the panel defaults to an area constraint. Proposed: `ddg_tools.surface_intersects(positions, faces)`.
* **Performance (WP17, `NOTES-WP17.md`, `benches/`):** `BlockTree` (BVH-ordered near blocks,
  parallel symmetric build, AVX near-field rows; products bitwise equal), batched GMRES (Schur
  columns + direction in one pass), 3-column `L′` solves, parallel GMRES vector operations. Flow
  defaults changed with measurements: GMRES tolerance 1e-3 (was 1e-4; 11 scenes × 200 steps, every
  step accepted, energies within trajectory noise, ≈ 20 % fewer iterations; still below the
  compressed operator's own error) and restart 60; `reference_compat` keeps 1e-4 / 30. Fig. 13 at 4
  threads: 12 of 14 models at or below the paper; above: #3 broken figure-8 0.50 s vs ≈ 0.37 s, #14
  figure-8 tunnel 3.15 s vs ≈ 2.7 s (remeshed to 53.5k faces, 59 GMRES iterations; the reference on the
  same mesh and machine 4.0–5.2 s) — **the 2.5 s at 45k target is not met** (machine load 6–16 from
  outside jobs during the timings). 1 thread vs the reference: faster on 5 of 8 scenes (walnut 140 vs
  259 ms, lung 482 vs 506, bunny 3522 vs 4205, reconstruction 7638 vs 9764, blocky Klein 237 vs 403),
  slower on EDT 1.09×, handcuffs 1.18×, Hopf 1.45× (ours 30/27/24 GMRES iterations vs the reference's
  21/13/9, which stop on the preconditioned residual). Thread scaling (40k torus): 3.32 / 2.08 / 1.44
  / 1.14 / 0.97 / 0.93 / 1.03 s at 1–64 threads; past 16 threads the back-substitutions, remeshing and
  Barnes–Hut gradient passes are sequential. The obstacle term in the preconditioner does not help
  (option, off).
* **Flow fidelity (`NOTES-FIDELITY.md`):** WP15's "torus drift" and symmetry table came from a bug in
  the example's `eig3` (fixed); a torus of revolution stays one within 0.6 % (reference 0.78 %);
  ablation of the symmetry error per layer: exact + dense 7e-15, hierarchical 2e-6, Barnes–Hut
  θ = 0.25 / 0.5 2e-4 / 1.6e-3, remeshing 4e-3. Remeshing every step sets an energy noise floor
  (tangential smoothing undoes ≈ the step's decrease near a minimum; polishing without remeshing
  lowers `E` by 2–6 %). Genus-g: 800 steps are too few (the linear arrangement is a slow saddle);
  g = 3 converges in ≈ 2080 steps to a near-tetrahedral shape (E 39.91; reference at 1000 steps
  39.54), g = 2 to E 38.59 (reference 38.40, same shape). Handcuffs unlink with and without CCD, in
  both programs, without the plane attractor; a strong plane (our weight 5 ≈ the reference's 150; the
  reference's plane is not area-weighted) keeps them linked — use weight ≤ 0.5. Eversion: the setup
  needs a gradually raised volume weight ≈ ℰ₀/V₀² and a symmetry search over principal axes; not
  reproduced yet. Bug fixed: `ccd::intersecting_pairs` flagged coplanar neighbours (orientation
  noise ±1e-20) — embedded meshes got a non-empty CCD band.
* **Solver robustness (in progress, `NOTES-SOLVERFIX.md`):** preconditioner scale `L_{2−s}/(w_TP ·
  4 I(σ) I(2−s))` (pins: GMRES stalled at 0.965 → converges), coarse correction on per-component
  translations (cobblestones 2800 → 134 iterations per step), remeshing growth cap and projection
  admissibility (remesh explosion and slab breach fixed); open: Klein bottle remeshing inside the
  CCD band at step ≈ 184, needle inputs of the subcritical spindle, envelopes (obstacle not in the
  preconditioner), punctured torus diverges after 38–80 steps (WP17), open surfaces without
  remeshing need 225–678 GMRES iterations.
* Fig. 16 (Ishii et al. 2012 table): no arXiv version exists; the user supplied the PDF
  (`/var/data/discrete-geometry/ishii2012.pdf`, not in the repository); transcription in progress. WP12 read "figure-8 / trefoil tunnel" as knotted
  tori (tubes around the knot); recorded as an interpretation.

## 17. The user's decisions (2026-11-01)

* Core in full (energy, `Hs` metric, Barnes–Hut, hierarchical matrices + GMRES with the
  two-Laplacian preconditioner, Schur and positional constraints, remeshing, every Sec. 7
  constraint and penalty).
* Applications: all of Sec. 9 (mathematical visualization, shape modeling, generative growth,
  subcritical `p < 4`).
* Evaluation: all six comparison methods, Fig. 5 / 10 / 12 / 13.
* Collisions: CCD in the line search, on by default only if measured cheap (D20).
* Blender: as for curves (N steps + live), pins from a vertex group untouched by remeshing,
  obstacles from objects, implicit obstacles/attractors from empties; remeshing replaces the mesh,
  no attribute transfer.
* The MIT reference may be read and ported, with its notice and credit to all four authors; its
  submodules and MKL are fetched for the black-box build (as for curves).
