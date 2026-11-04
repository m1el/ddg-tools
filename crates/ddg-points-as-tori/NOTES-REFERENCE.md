# Points as Tori: notes on the reference implementation

Survey of 2026-11-04. This file records what an implementer needs and the paper does not say exactly.

**Paper.** Nicole Feng, Ioannis Gkioulekas, Keenan Crane. *Points as Tori: Fast Pointwise Signed
Distance for Point Clouds.* ACM TOG 45(4), Article 53, 2026 (`../PointsAsTori/PointsAsTori.pdf`).

**Reference.** `../PointsAsTori/points-as-tori` at `253bf9f` (2026-07-29), MIT, "Copyright (c) 2026
Nicole Feng". It may be read and ported; the notice is kept in `LICENSES/points-as-tori-MIT.txt`.
Submodules (`deps/{nanobind,libigl,fcpw,nanoflann}`) are empty and are not fetched here.

**Citation forms.**
- `ref: <file>:<line>` points into the reference.
- `NW` = `src/pointsastori/network.py`, `TD` = `src/cpp/torus_distance.{h,cpp}`,
  `S3` = `src/cpp/shape_3d.cpp`, `S2` = `src/cpp/shape_2d.{h,cpp}`, `T` = `training/train.py`,
  `MS` = `training/mesh_sampling_jax.py`, `MD` = `training/mesh_data.py`,
  `SH` = `demo/shaders_3d.py`, `DM` = `demo/demo_3d.py`.
- Equation numbers are the paper's. "Verified" means checked by a script in this survey;
  the scripts are listed in §20.

## 0. Summary of surprises

1. **The 2D "circles" variant (Fig. 10) does not exist in the reference.** `shape_2d.*` is a
   library of analytic Quilez 2D SDFs, used to make the 3D test shapes by extrusion and
   revolution. Fig. 10's caption says "In 2D, our method *would* use circles". There are no 2D
   weights and no 2D network (§13).
2. **Normal orientation (Fig. 27) is not in the reference** (§15).
3. **The released evaluation is not the one in Sec. 4.2.**
   - `infer.py` evaluates with the 32 nearest points and a *per-query* `λ = 64/σ_x`, where `σ_x`
     is half the largest of those 32 distances.
   - The paper's scheme (`λ = 10³/D` constant, sum over `R_eval = 2C/λ`, 32-NN fallback) appears
     only as dead code: `compute_optimal_radius` and `set_radius_evaluation` are never called.
     Its `λ` rescaling (`λ *= R_eval/√3`) is dimensionally inconsistent (§7.3).
4. **Architecture.** The pretrained network is 8 layers, 8 heads, MLP 512, `embed_dim` 128, k = 64.
   This matches the paper and not the constructor defaults (4/4/256).
   - There is **no attention mask**.
   - GELU is the **tanh approximation**.
   - LayerNorm uses **ε = 1e-6** and the "fast" variance `E[x²] − E[x]²`.
   - The query is pre-scaled by `1/√16`.
5. **Two different pretrained weight sets.** The repository's `models/FundamentalFormPredictor.pkl`
   and the training data's `FundamentalFormPredictor_pretrained.pkl` have the same configuration,
   but all 136 arrays differ (max |Δ| ≈ 0.18). The repository file is the default for inference.
6. **Training differs from Sec. 5.1 in many places** (§9.7): an unmentioned lr = 1e-3 pretraining
   phase, "cosine decay" that is a no-op, query box sizes, and more.
   - The network runs on all 2048 neighborhoods of a cloud at every step. The "batch of 1024" is
     a batch of *query* neighborhoods.
   - Padded clouds give NaN gradients, which are zeroed, so 3,588 of the 9,998 simulated shapes
     probably never contribute.
   - The pretraining blob and the validation files have **normals of length 2**.
7. **Precision.**
   - The network and all of training are f32: `jax_enable_x64 = False`, `NW:23`.
   - Fitting and evaluation in C++ are f64. The demo shader is f32 GLSL.
   - f32 features are computed as `positions/σ − p₀/σ` on *absolute* coordinates. This loses
     precision for clouds far from the origin.
   - Fitted tori of nearly flat regions have `R, r ~ 1/|κ|` of up to `1e10`, so the torus SDF
     `‖d‖ − |r|` cancels catastrophically in f32 (§18).
8. **No input validation** anywhere: NaN tori, zero normals, duplicate points, fewer than k
   points, and query = data point (`r̂ = 0/0`) all give NaN or garbage silently (§17).

## 1. Inference data flow (`ref: src/pointsastori/infer.py`)

`PointsAsTori(points, normals)` runs these steps:

1. `TorusDistanceField(points, normals)` (C++, f64) copies the cloud and builds a nanoflann kd-tree
   (`max_leaf = 20`). It also sets an SDF box of twice the cloud's bounding box, used only for
   sampling (`TD.h:155-183`).
2. Loads the model (`load_saved_model`, `NW:1174-1208`). This unpickles the pkl without
   restriction; we do not.
3. `precompute_coefficients_in_chunks(points, normals, k = 64, chunk_size = 50000)`
   (`NW:1122-1172`):
   - `get_neighbors` (C++ nanoflann, §4).
   - Gathers `points[idx]`, `normals[idx]` and converts them to f32 jnp arrays.
   - Runs the network forward in chunks of 50k neighborhoods.
4. `fit_tori_from_forms(points, normals, coeffs)` (C++, f64, `TD.cpp:779-837`, §6). The f32
   coefficients are promoted to f64.
5. `signed_distance(queries, accelerate = True)` calls `set_k_evaluation(32)` and then
   `evaluate_distance` (§7). `accelerate = False` sums over all points (`k = -1`).
6. `sdf_gradient` returns the analytic gradient (§7.4).
7. `save_tori` and `_load_tori` pickle the tuple `(centers (n,3), axes (n,2), R (n,), r (n,))`.
   The sign of the torus lives in the sign of `r`, and the axis is stored as spherical angles
   `(θ, φ)`.

The point cloud is used as given: no centering, scaling or normal normalization. `read_point_cloud`
reads PLY through open3d (normals not normalized). For OBJ it normalizes `vn` (`shape_3d.py:42-72`).

## 2. Network (`NW:910-1110`)

`FundamentalFormPredictor(embed_dim = 128, n_freq = 0, n_layers = 8, n_heads = 8, mlp_dim = 512,
dropout_rate = 0.0)`. The constructor defaults are 4/4/256; the pkl stores 8/8/512, which matches
the paper's Sec. 4.3. Flax NNX 0.12.2, JAX 0.8.2 (`pyproject.toml`).

### 2.1 Input

`positions, normals: [B, k, 3]`. Row 0 is the center point `p_i`; rows 1..k-1 are its neighbors
in ascending distance (§4). The neighborhood **includes `p_i`**.

### 2.2 Features (`InputFeatures.__call__`, `NW:983-1028`), all f32

| quantity | formula |
|---|---|
| unit normals | `n̂_j = n_j / max(‖n_j‖, 1e-10)` |
| center normal | `c = n̂_0` |
| `σ_i` | `jnp.median_j ‖p_j − p_0‖` over all k = 64 entries, **including `j = 0` (distance 0)**: the mean of the 32nd and 33rd smallest of the 64 values (indices 31 and 32 when sorted ascending) |
| scaled positions | `p'_j = p_j / σ_i`, **absolute coordinates divided first**, then `d_j = p'_j − p'_0` (`NW:1006-1013`) |
| local basis | `u = (0,1,0)` if `|c_x| > 0.9` else `(1,0,0)`; rotation 0, so `u_rot = u` (Rodrigues with θ = 0 is the identity up to round-off: `cos 0 = 1`, `sin 0 = 0` exactly); `ŝ = (c × u)/max(‖c × u‖, 1e-10)`; `t̂ = c × ŝ` (`NW:35-80`) |
| features | `[d_j·ŝ, d_j·t̂, d_j·c, n̂_j·ŝ, n̂_j·t̂, n̂_j·c]` in this order, ∈ R⁶ |

- `n_freq = 0`. The code path for `n_freq > 0` changes `input_dim` to `3(2n+1)+3` but never
  computes an encoding (`NW:970-973`), so it would crash on a shape mismatch. Dead code.
- The C++ `orthonormal_basis` (`S2.cpp:20-33`) is the same basis: `n/‖n‖` with no ε,
  `AngleAxis(0)`, then `s/‖s‖`. It is used in torus fitting, so network features and fitting share
  the frame up to round-off. That only matters if `‖n‖` differs between the two; both normalize.

### 2.3 Layers

```
x = LayerNorm₁₂₈(Linear₆→₁₂₈(features))                 # input_encoder.mlp.layers.{0,1}
repeat 8×:                                              # blocks.{0..7}, pre-LN
    x = x + MHA(LN₁(x))                                 # no dropout effect (rate 0)
    x = x + Linear₅₁₂→₁₂₈(GELU_tanh(Linear₁₂₈→₅₁₂(LN₂(x))))   # ffn.layers.{0,3}; 1, 2, 4 are gelu/dropout
x = LayerNorm(x)                                        # final_norm
y = Linear₁₂₈→₆(x[:, 0, :])                             # mlp.layers.0, center token only
coeffs = (y0·σ, y1, y2, y3/σ, y4/σ, y5/σ)               # (a00, a01, a10, a11, a02, a20)
```

**Linear.** `y = x @ W + b`, with `W: [in, out]` (Flax layout).

**LayerNorm** (Flax `nnx.LayerNorm`, verified from the flax 0.12.2 source):
- `μ = mean(x)`, `var = max(0, mean(x²) − μ²)` ("fast variance").
- `y = (x − μ)·rsqrt(var + 1e-6)·scale + bias`, with `ε = 1e-6`.
- Statistics are computed in at least f32.

**GELU.** `nnx.gelu = jax.nn.gelu(approximate = True)`, i.e.
`0.5x(1 + tanh(√(2/π)(x + 0.044715x³)))` (verified).

**MHA** (`nnx.MultiHeadAttention(num_heads = 8, in_features = 128, qkv_features = 128,
out_features = 128)`):
- `q, k, v = LinearGeneral(x)` with kernel `[128, 8, 16]` and bias `[8, 16]`, so the head
  dimension is 16 and the biases are on. `normalize_qk = False`.
- `q ← q/√16`, applied to the query before the dot product, as in Flax `dot_product_attention`.
- `A = softmax_k(q·kᵀ)` per head, with **no mask**, all 64 tokens.
- `o = Σ A v`, then `out = LinearGeneral` with kernel `[8, 16, 128]` and bias `[128]`.
- Softmax is in f32.
- Padded neighbors (training only, §9.2) are not masked. They are real tokens whose features come
  from the zero row.

**Output scaling** (`NW:1097-1108`) matches the paper. a00 is a length; a01 and a10 are slopes;
a11, a02 and a20 are curvatures (1/length).

**Init** (`NW:1058-1067`).
- The final kernel is zero and the final bias is `(0, 0, 0, 0, −0.5, −0.5)`.
- So the initial `Q(s,t) = −(s² + t²)/(2σ)`. This is the sphere of radius σ tangent at `p_i`,
  inside, with the normal pointing out (§6.3).
- All other parameters use Flax defaults:
  - Linear and MHA kernels: lecun_normal (variance_scaling 1, fan_in, truncated normal).
  - Biases: 0.
  - LayerNorm: scale 1, bias 0.

### 2.4 Parameter names and shapes (both pkls, f32)

| name | shape |
|---|---|
| `input_encoder.mlp.layers.0.{kernel,bias}` | `[6,128]`, `[128]` |
| `input_encoder.mlp.layers.1.{scale,bias}` | `[128]`, `[128]` |
| `blocks.{L}.norm1.{scale,bias}`, `blocks.{L}.norm2.{scale,bias}` | `[128]` each |
| `blocks.{L}.attention.{query,key,value}.kernel` | `[128,8,16]` |
| `blocks.{L}.attention.{query,key,value}.bias` | `[8,16]` |
| `blocks.{L}.attention.out.kernel` / `.bias` | `[8,16,128]` / `[128]` |
| `blocks.{L}.ffn.layers.0.{kernel,bias}` | `[128,512]`, `[512]` |
| `blocks.{L}.ffn.layers.3.{kernel,bias}` | `[512,128]`, `[128]` |
| `final_norm.{scale,bias}` | `[128]` |
| `mlp.layers.0.{kernel,bias}` | `[128,6]`, `[6]` |

- `L = 0..7`. There are 136 arrays and 1,588,358 parameters (6,353,432 bytes).
- Top-level pkl keys: `params` (an `nnx.State` of `Param`s), `embed_dim = 128`, `n_freq = 0`,
  `n_heads = 8`, `n_layers = 8`, `mlp_dim = 512`, `k_neighbors = 64`, `dropout_rate = 0.0`.
- Globals in the pickle: `flax.nnx.statelib.State`, `flax.nnx.variablelib.Param`,
  `flax.nnx.tracers.TraceState`, `jax._src.array._reconstruct_array`,
  `numpy._core.multiarray._reconstruct`, `numpy.ndarray`, `numpy.dtype`.
- The pkl was read with a restricted `Unpickler.find_class` that maps exactly these to stubs
  (§20); the `TraceState` and metadata are discarded.
- The trained final bias is `(0.014, 0.006, 0.011, 0.003, −0.482, −0.492)`, close to init. Its
  kernel max is 0.164.

## 3. Weights files

- `src/pointsastori/models/FundamentalFormPredictor.pkl`: 6,366,869 bytes. The default for
  inference (`infer.py:11`).
- `PAT_training_data/FundamentalFormPredictor_pretrained.pkl`: same size and config, but a
  **different** set of weights.
  - All 136 arrays differ; for example, `blocks.0.attention.key.kernel` has max |Δ| = 0.18.
  - The name suggests a checkpoint of another phase, perhaps the blob pretraining or a
    different run.
  - Nothing says which phase it is. Convert both and record which is which.

## 4. Neighbors (`TD.h:51-117`, nanoflann)

- Uses `KDTreeEigenMatrixAdaptor<Matrix<double,3,Dynamic>, 3, metric_L2_Simple>` with
  `max_leaf = 20` and an exact kNN (`KNNResultSet`). Results are sorted ascending by squared
  distance.
- **The query point is its own first neighbor** (distance 0). Verified on all valid training rows
  (§9.2).
- Ties, including duplicate points, are ordered by kd-tree traversal, which is unspecified.
  - With duplicates, `p_0` may not be the query point.
  - The JAX feature code assumes row 0 is the center. With a duplicate first, the center's own
    normal could be replaced by the duplicate's.
- `k ≥ n` gives k = n; the matrix is then `(n, n)` and the network sees fewer tokens. For training
  data, columns are padded with −1 (§9.2).
- With outliers (`outliers` argument): the tree is built on the inliers, and **outlier columns
  stay −1**. In Python, `points[-1]` silently picks the last point. This is a bug; nobody calls it
  with outliers.
- **Evaluation-time kNN** (`TD.cpp:76-151`) reuses the same tree for arbitrary queries.
  - k = 32 (`infer.py:12`).
  - The radius search (`radiusSearch` on squared radius, unsorted) is unused at inference.
- `get_average_neighbor_distance` (`TD.h:119-150`) has a **data race**: `total_distance +=` inside
  `omp parallel for` with no reduction. It only feeds the unused `compute_optimal_radius`.

## 5. Coordinates and symbols

- `Q*_i(s,t) = p_i + s ŝ_i + t t̂_i + Q_i(s,t) n_i`, with
  `Q_i = a00 + a10 s + a01 t + a11 st + a20 s² + a02 t²` (Eq. 18; `a_{n,m}` multiplies `sⁿ tᵐ`).
- Code order: `coeffs = (a00, a01, a10, a11, a02, a20)`.
- **`n_i` is the raw input normal** in fitting (`TD.cpp:652-653`, `NW:212-213`). It is not
  normalized; the basis is built from the normalized normal. See §17, item 4.

## 6. Torus fitting

### 6.1 As coded (C++, `get_torus_parameters_python` + `fit_torus_at_point`, `TD.cpp:647-807`; JAX `torus_precompute_masked`, `NW:183-336`, identical math)

```
A    = √(1 + a01² + a10²)
v1   = a20 a01² − a02 a10² − a02 + a20
v2   = a11 (a10² + 1) − 2 a01 a10 a20
term = a02 (a10² + 1) − a01 a10 a11 + a20 (a01² + 1)
disc = √max((1 + a01² + a10²)(a11² − 4 a02 a20) + term², 1e-10)          # NAN_EPSILON clamp
ev₊  = (sgn⁺(v1 + disc)·max(|v1 + disc|, 1e-10), v2)
ev₋  = (v1 − disc, sgn⁺(v2)·max(|v2|, 1e-10))
ev±  /= ‖ev±‖
I    = [[1 + a10², a10 a01], [a10 a01, 1 + a01²]]
II   = [[2 a20, a11], [a11, 2 a02]] / A
κ±   = (ev±ᵀ II ev±) / clampsign(ev±ᵀ I ev±, 1e-10)                        # Rayleigh quotients
v±   = ev±_x (ŝ + a10 n) + ev±_y (t̂ + a01 n)                               # 3D directions, unnormalized
n*   = (n − a10 ŝ − a01 t̂) / A                                              # Eq. 21
if |κ₊| < |κ₋|: κmin = κ₊, κmax = κ₋, vmin = v₊, else κmin = κ₋, κmax = κ₊, vmin = v₋
r    = 1/|κmax|;  R = 1/|κmin| − sgn(κmax κmin) r                           # sgn(0) = 0 here
s    = +1 if (κmax < 0 ∧ κmin < 0) ∨ (κmax < 0 ∧ κmin > 0) else −1          # Eq. 22
c    = p + a00 n − (s/|κmin|) n*
u    = (n* × vmin)/max(‖n* × vmin‖, 1e-10)
axis = (atan2(u_y, u_x), acos(u_z))                                        # stored; C++ does not clamp u_z
torus = (c, axis, R, s·r)                                                  # sign stored in the minor radius
```

**Sign helpers.**
- `sgn⁺(x) = +1` for `x ≥ 0`, else −1 (C++ `sign`, `S2.h:25`; NaN gives −1).
- JAX uses `jnp.sign`, which is 0 at 0. This differs only at exact zeros.
- `sgn` in `R` is C++ `sgn` (`S2.h:21`, can be 0), matching `jnp.sign`.

**Clamps.**
- JAX additionally clamps `|κmin|, |κmax| ≥ 1e-10` before inverting (`NW:306-311`).
- The C++ inference version does not (`TD.cpp:797-801`): exactly zero curvature gives
  `r = ∞` or `R = ∞`, which leads to NaN tori.
- The C++ has no NaN filter at evaluation; JAX masks non-finite tori in training (`NW:394-401`).

**`acos` domain.** The JAX `cartesian_to_rotation` clamps `u_z` to [−1, 1] and guards `x ≈ 0`
(`NW:88-112`). The C++ `cartesian_to_rotation` (`TD.h:36-40`) does not; `|u_z| > 1` by round-off
gives NaN.

**Eigenvectors.**
- `(v1 ± disc, v2)` are generalized eigenvectors of `(II, I)`, and the Rayleigh quotients equal
  the paper's `κ± = H ± √(H² − K)` (Eq. 19).
- Verified on 10⁵ random coefficient vectors: eigen-residual `< 6e-12` and mismatch `< 6e-14`
  relative (§20).
- This is a different (equivalent) formula from Eq. 20's `w± = (κ± a10 a01 A − a11,
  2 a20 − κ±(1 + a10²) A)`. Both come from one row of `(II − κ I) w = 0`.
- The reference's other copies, `get_torus_parameters` and `get_polynomial_curvatures`, clamp disc
  at 0 instead of 1e-10 and have no eigenvector clamps. They are used only by
  `principal_curvatures`/`principal_directions`, which are exported but unused.

### 6.2 Paper vs code (Sec. 4.1)

| | paper | code |
|---|---|---|
| curvatures | `H ± √(H² − K)`, Eq. 19 | Rayleigh quotients of the code's eigenvectors (equal) |
| `v̂min` | normalized `v±` | unnormalized; only `n* × vmin` is normalized (same axis) |
| clamps | none | 1e-10 clamps on `disc`, `ev₊[0]`, `ev₋[1]`, `ev·I·ev`, and (JAX only) `|κ|` |
| `R` | `1/|κmin| − sign(κ₊κ₋) r` | same, with `sgn(0) = 0` |
| sign | Eq. 22 | same expression; equals "`κmax < 0` and `κmin ≠ 0`" |

### 6.3 Degenerate cases (what the code produces)

- **Sphere / umbilic** (`a20 = a02`, `a11 = a10 = a01 = 0`).
  - `v1 = v2 = 0` and the true disc is 0, clamped to `1e-5`.
  - `ev₊ = (1, 0)`. `ev₋ = (−1, 0)` in JAX; `(−1, 1e-5)` normalized in C++.
  - Both curvatures equal `2a20`, so `κmin = κ₋` (strict `<` fails), giving `R = 0`, `r = 1/|κ|`,
    `s = +1` for `a20 < 0`.
  - The axis is arbitrary but irrelevant when `R = 0`. The initial network output is exactly this
    case with `κ = −1/σ`: center `p − σn`, radius σ.
- **Plane** (all second-order coefficients 0).
  - JAX: `|κ| → 1e-10`, so `r = R = 1e10`, `s = −1`, `c = p + 1e10 n*`.
  - C++: `r = ∞` gives `R = NaN` (`∞ − 0·∞`) and a NaN center. The network never outputs exact
    zeros, so in practice `R, r` can be `~10³–10⁸` for flat regions.
  - The SDF `‖d‖ − |r|` then cancels catastrophically (§18).
- **Cylinder** (one κ = 0).
  - `κmin → 0`, so `R → ∞` and `c` goes to infinity along `n*`. The axis `n* × vmin` is the
    cylinder's generator direction.
  - With JAX clamps, `R = 1e10 − sgn(0)r = 1e10`. Same cancellation problem.
- **Saddle** (`κ₊κ₋ < 0`).
  - `R = 1/|κmin| + r`, `s = +1` iff `κmax < 0`.
  - The torus's inner equator touches `Q*(0,0)`; this is the paper's "κmax < 0, κmin > 0" case
    (Fig. 14).
- **`R < r`.** "We take `R` always greater than `r`" (paper) holds by construction for
  `κ₊κ₋ > 0` (`R = 1/|κmin| − 1/|κmax| ≥ 0`). `R ≥ r` does not hold in general: a self-intersecting
  spindle torus results when `|κmin| > |κmax|/2`. The SDF formula is still the "circle twice"
  formula.
- **Zero or non-unit normal.**
  - Zero: basis `ŝ = 0/1e-10 = 0` (JAX) or NaN (C++), so garbage.
  - Non-unit: features use the normalized normal, fitting uses the raw one (§17, item 4).

## 7. Evaluation (C++, `TD.cpp:153-534`)

### 7.1 Per-torus SDF (Eq. 23–24, `TD.cpp:153-216`)

```
a = (cos θ sin φ, sin θ sin φ, cos φ)                   # from the stored angles (TD.h:29-34)
x = q − c;  v = ‖x × a‖;  d = (v − R, x·a)
g(q) = sgn⁺(r)·(‖d‖ − |r|)
∇g   = sgn⁺(r)·(w₁ G₁ + w₂ a),  w = d/‖d‖,  G₁ = ∇v = (x − (x·a)a)/v
```

- `G₁` is written componentwise in `TD.cpp:199-207` as the negated, expanded form of
  `−(I − aaᵀ)(c − q)/v`.
- There are **no guards**: `v = 0` (q on the axis) and `‖d‖ = 0` (q on the core circle) give NaN
  in C++. JAX uses `max(·, 1e-10)` (`NW:448-475`).
- Storing the axis as angles and recomputing `a` costs `sin`/`cos` per evaluation, and its round-off
  is not negligible in f32. We store `a` as a unit vector (D-list in SPEC).

### 7.2 Blend (Eq. 25–26 with the released neighbor rule)

For a query `q` and a neighbor set `𝒩(q)`:

```
σ_q = ½ max_{i∈𝒩} ‖q − p_i‖              # "shift", compute_shift (TD.cpp:17-56)
λ_q = 64/σ_q · lambda_scale               # compute_lambda (TD.h:42-44); lambda_scale = 1 by default
w_i = exp(−λ_q(‖q − p_i‖ − σ_q)) [· A_i if use_areas]
φ(q) = Σ w_i g_i(q) / Σ w_i − isovalue
```

**Neighbor sets.**
- Default `accelerate = True`: the 32 nearest data points to `q`.
- `accelerate = False`: all points. Then `σ_q` is half the distance to the farthest point of P,
  which is Eq. 26 exactly.
- Radius mode: unused (§7.3).

**Weight range.** With 32-NN, `λσ = 64` and `d_i ∈ [0, 2σ]`, so the exponent is in [−64, 64]. This
is the "C = 64, below 87 for f32" argument. For all points, `d_i ≤ 2σ` also holds by definition
of `σ` as half the maximum.

**Area weights.** `use_areas = false` by default. The paper says it omits area weights, and the
demo shader uses them in its brute-force mode (§14).

**Outliers.** Points with `is_outlier` are skipped in the sum, but still count in `σ_q`.

**Isovalue.** It is subtracted at the end, so offsets are `φ − c`, with no refit.

### 7.3 Radius mode (Sec. 4.2 as printed) — dead code

- `set_radius_evaluation(r)` switches to `radiusSearch(r)`. If no point is found, it falls back to
  kNN with `k_evaluation` or 10.
- It then **multiplies `λ` by `r/√3`** "for radius-based truncation" (`TD.cpp:274-280`; the comment
  says "I can't remember why I did this"). That makes `λ` dimensionless, `≈ 128/√3` when `σ ≈ r/2`.
- `compute_optimal_radius(points, 64)` (`shape_3d.py:1193-1196`) computes `D` = the average
  64-NN distance (racy, §4), `λ = 10³/D` and returns `128/λ`. This is the paper's `R_eval = 2C/λ`.
- Nothing calls either function. The paper's constant `λ = 10³/D` is never used at evaluation: the
  per-query `λ = 64/σ_q` is used in all modes.

### 7.4 Gradient and Laplacian (`TD.cpp:319-408`; App. D)

```
W = Σ w_i;  r̂_i = (q − p_i)/‖q − p_i‖
A = Σ w_i (∇g_i − λ r̂_i g_i)/W;  B = λ Σ w_i r̂_i / W
∇φ = A + φ B
```

- This treats `σ_q` and `λ_q` as constants. They are piecewise smooth functions of `q`, and the
  printed gradient ignores their derivative. That is exact for `λ` and `σ` fixed. With
  `λ_q = 64/σ_q` varying, the true gradient has an extra term `(∂λ/∂q)·Σ w_i(σ − d_i)(g_i − φ)/W`
  (the shift cancels).
- So **`∇φ` is not exactly the gradient of `φ`** where the farthest neighbor changes smoothly.
  The error is small but measurable, and it matters for finite-difference tests (SPEC decision).
- `r̂` with `q = p_i` exactly is `0/0`, so NaN in C++ (JAX guards). This happens when evaluating at
  the data points (normal orientation, morphology).
- **Laplacian**: the App. D divergence formula, with `Δg_i` replaced by a placeholder
  `s (d − 1)/‖q − c_i‖` ("TODO", `TD.cpp:392`). The README lists it as unfinished. It is wrong for
  every torus; not ported.

### 7.5 Batch evaluation

- Queries are parallel with OpenMP (`schedule(static)`, or `dynamic` for radius mode).
- Each query is independent, so results are deterministic.
- The kNN is computed for all queries first (`compute_nearest_points`), stored as a `k × n_queries`
  matrix, and then evaluated.

### 7.6 Other evaluators in `TD.cpp` (baselines)

| function | formula | notes |
|---|---|---|
| `signed_log_sum_exp` (`TD.cpp:987-1060`) | SHC, Eq. 28: `w = Σ A_i (λρ+1)⟨q−p_i, n_i⟩/(2πρ³) e^{−λ(ρ−σ)}`, `φ = sgn(w)(−log|w|/λ + σ)` | needs `areas`; same `σ`, `λ = 64/σ` per query and kNN rule as PAT |
| `signed_log_sum_exp_gradient` | closed form | the `term2` coefficient `6πρ…/(4π²ρ⁶)` should be checked against `∂(ρ⁻³)`: `6π r r_vec/(4π²r⁶)` = `3 r_vec/(2π r⁵)`, correct |
| `self_normalized_signed_log_sum_exp` (`TD.cpp:1118-1138`) | SSPD, Eq. 29 | needs `areas` |
| `log_sum_exp` (`TD.cpp:1062-1077`) | Madan–Levin unsigned LSE with `λ = lambda_scale`, `shift = 64/λ` | all points |
| `winding_number` | `Σ A_i ⟨n_i, q−p_i⟩/(4πρ³)` | brute force |
| `regularized_winding_number` | Chen et al. 2024: `S(t) = erf(t) − 2t/√π e^{−t²}`, `t = ρ/ε` | brute force |
| `fast_winding_number_precompute` / `fast_winding_number` | libigl octree, `expansion_order = 2`, `β = 2` | **bug**: the precompute computes point areas into a local `A` (`TD.cpp:1272-1279`) and never stores them; the query passes the member `areas`, which is empty unless `get_point_areas()` ran first |
| `compute_point_areas` | `igl::octree` + `igl::knn(k = 20)` + `igl::copyleft::cgal::point_areas` (Barill et al. 2018, Sec. 5) | the FWN path uses k = 21 ("used in CGAL fast_winding_number"); CGAL is GPL, so we implement Barill's construction from the paper |

## 8. Ground truth and meshes (`S3.cpp`)

- `Mesh3D` uses FCPW for unsigned distance and closest point (`S3.cpp:712-760`). The sign comes
  from the generalized winding number, `evaluate_gwn`, brute force (`S3.cpp:854-893`), and a
  libigl FWN variant exists.
- `remove_degenerate_faces` drops faces with area `< 1e-10` (`S3.cpp:788-815`).
- `center_and_scale` subtracts the **vertex mean** and divides by `max|coord|`
  (`S3.cpp:54-61`), so the shape is in [−1, 1]³ and touches it.
  - The paper's evaluation says "centroids at the origin, scaled to lie within [−1, 1]³".
  - Training uses the C++ version; the paper's Thingi10K scaling is unspecified.
- In training, ground truth is recomputed in JAX (`MS:16-151, 395-424`):
  - exact point–triangle distance,
  - sign from the GWN sum `Σ 2 atan2(det[a,b,c], abc + (a·b)c + (b·c)a + (c·a)b)`,
  - inside iff `gwn > π` (that is, winding number > 1/4).

## 9. Training (`T`, `MS`, `MD`)

### 9.1 Loss

`SDF_loss`, `T:59-122`; `torus_distance_gradient_masked`, `NW:367-518`.

**Per step:**
1. The network runs on **all** 2048 neighborhoods of one cloud (`T:246`).
2. All 2048 tori are fitted (`NW:183-336`, f32, with clamps, §6.1).
3. The loss is evaluated on a contiguous slice of `batch_size` neighborhoods' queries
   (`dynamic_slice`, `T:510-514`).
4. `batch_neighbors` and `batch_key` are computed but unused (`T:507, 511`).

**Prediction** is Eq. 25 over the whole cloud's valid tori:
- `σ_q` = half the max distance among the **64** nearest valid points to `q`
  (`lax.top_k` over 2048, `NW:142-166`).
- `λ_q = 64/σ_q`.
- `w_i` is masked for invalid or non-finite tori.
- `φ = Σ w g/max(W, 1e-10)`.
- `∇φ` is the closed form of §7.4 with 1e-10 guards.
- The weights depend only on `q` and the `p_i`, not on parameters. Parameters reach the loss through
  `g_i` and `∇g_i`.
- Training therefore uses **64-NN** for `σ` (inference uses 32), and sums over all tori (inference
  sums over 32).

**Loss:**
- `L = mean|φ − φ_true| + mean|1 − ‖∇φ‖|`, both L1, weight 1.
- The means are over valid (row, query) pairs; NaN predictions are dropped from the mean.
- No regularizer or weight decay.
- NaN gradient entries are set to 0 after `value_and_grad` (`T:254-260`).
- The eikonal term differentiates through the closed-form `∇g`, so it needs second derivatives of
  the torus SDF with respect to (q, torus parameters).

**Queries** (`MS:609-681`):
- Regenerated per shape per epoch: 120 per neighborhood, 40 + 40 + 40, for all 2048 rows including
  padding.
- Per-row key: `split(shape_key, 2048)[i]`.
- Row box: the AABB of the row's 64 neighbors. A regularized box with each extent `≥ 0.1` is used
  for face selection (`MS:564-569`).

| set | count | how |
|---|---|---|
| zero set | 40 | faces whose AABB meets the regularized box, chosen ∝ area; barycentric `(1−√u, √u(1−v), √u v)`-style (`MS:286-305`); can leave the box |
| narrow band | 40 | same faces (independent key), offset along the face normal by `U(−0.2, 0.2)` (`MS:495-523`; `T:307`); absolute units |
| cube | 40 | uniform in a cube centered at the row AABB center, half-side `0.5·√3·max(max extent, 0.01)` (`MS:595-601`; `bbox_extension = 1`) |

### 9.2 Batching, masks, RNG

**Steps per shape.** `n_batches = ceil(2048/batch_size)`:
- 2 steps per shape at batch 1024; both steps share the queries.
- 64 steps per shape at batch 32 (pretraining).
- Steps over all-padding rows still step Adam with a zero loss.

**Arrays.**
- `masks (S, 2048)`: a valid prefix.
- `neighbors (S, 2048, 64) int32`:
  - Column 0 is the point itself (verified).
  - Padded rows are all −1.
  - For `n < 64`, columns are −1-padded.
  - JAX indexes −1 as row 2047, the zero padding row whenever `n < 2048`. So padded tokens are
    origin points with zero normals, and **there is no mask in attention**.

**Shuffle and seeds.**
- `np.random.seed(42)`, then a new permutation of all shapes of the phase each epoch (`T:447, 458`).
- `PRNGKey(42)`, `fold_in(key, epoch)`, `fold_in(·, global_idx)` (`T:292, 450, 477`).
- `nnx.Rngs(42)` for init (`T:293`).
- Validation keys: `fold_in(epoch_key, 10000 + i)`.

**Reproducibility.** File lists come from unsorted `glob` (except phase 0), so reproducibility
depends on the filesystem order.

### 9.3 Curriculum (`train()`, `T:705-802`)

Each phase starts a new model, loads the previous phase's pkl, and gets a fresh Adam.

| phase | data | shapes | epochs | lr | warmup | batch | steps |
|---|---|---|---|---|---|---|---|
| pretrain | `pretraining_blob.npz` | 1 | 64 | 1e-3 | 0 | 32 | 4,096 |
| 0 | `nice/blobs_train_mode0*` (sorted) | 2000 | 24 | 1e-5 | 100 | 1024 | 96,000 |
| 1 | blobs + ABC17 mode0, `zip` (truncates to 4 ABC17 files) | 4000 | 8 | 1e-5 | 1000 | 1024 | 64,000 |
| 2 | blobs + all ABC17 + ABC2 | 7459 | 8 | 1e-6 | 1000 | 1024 | 119,344 |
| 3 | blobs + all `ABC*_train_mode0_*` | 12459 | 5 | 1e-6 | 1000 | 1024 | 124,590 |
| 4 | everything nice + all simulated | 22457 | 2 | 1e-7 | 1000 | 1024 | 89,828 |

Total ≈ 497,858 steps. The paper reports 45 h on a 3090, about 0.33 s per step.

### 9.4 Optimizer and schedule (`T:178-216`)

- `optax.adamw(warmup_cosine_decay_schedule(0, lr, warmup, decay_steps = 1e8, end = 0.01 lr),
  weight_decay = 0.0)`, which is Adam with β = (0.9, 0.999), ε = 1e-8.
- `decay_steps = 1e8` with at most 125k steps per phase means **linear warmup, then a practically
  constant lr**. The paper's cosine decay is a no-op.
- No clipping (commented out), no EMA, no gradient accumulation.
- A checkpoint is saved every epoch to `models/FundamentalFormPredictor_{epoch}.pkl`; names
  collide across phases.

### 9.5 Precision

- Training is f32. Before commit `c9e2e37`, `MS` enabled x64 on import, so older runs may have
  computed in f64.
- The 45 h figure implies f32: the f64 rate of a 3090 would need ≥ 370 h.

### 9.6 Validation

- `ABC17_validation_3.npz` (10 shapes, at the relative path `datasets/…`, `T:383`).
- Fresh queries each epoch; L1 distance and eikonal per batch; logged only.
- `ABC0_validation_0.npz` has an older format (P = 1024, **k = 128**) and is unused.

### 9.7 Paper vs code (training)

| item | paper (Sec. 4.3, 5.1) | code |
|---|---|---|
| `σ_x` in the loss | max over all of P (Eq. 26) | max over the 64 nearest |
| cube queries | side = largest extent of the neighborhood box | side = √3 × largest extent (floor 0.01) |
| narrow band | width 0.2 | ±0.2 (width 0.4), only from faces near the regularized box |
| FPS oversampling | 2048 · 8 | 2048 · 4 (`S3.cpp:163-203`) |
| scan sensor | side √3 | 2√3 (`sensor_size 1.0 · 2√3`) |
| learning rates | 1e-5 first phase, 1e-6, last 1e-7, warmup 1000 | as table above: phases 0 **and 1** at 1e-5, warmup 100 in phase 0, plus a pretraining phase at 1e-3 |
| schedule | linear warmup + cosine decay | warmup + constant (decay over 1e8 steps) |
| batch | 1024 neighborhoods | 1024 query rows; the network runs on all 2048 |
| epochs | 24, 8, 8, 5, 2 | same (plus pretraining 64) |
| ABC chunks | 2, 17, 55, 58 | same (files ABC2/17/55/58) |
| shapes per phase | 2000 blobs; +2000 ABC; +3500 ABC; all direct; all | 2000; 4000; 7459; 12459; 22457 |

## 10. Training data (`PAT_training_data/PAT_training_data/`, 5.7 GB)

**Members of every npz** (`MD:986-1018`):

| member | dtype | shape |
|---|---|---|
| `positions` | f64 | (S, P, 3) |
| `normals` | f64 | (S, P, 3) |
| `neighbors` | i32 | (S, P, k) |
| `masks` | bool | (S, P) |
| `mesh_vertices` | f64 | (S, V, 3) |
| `mesh_faces` | i32 | (S, F, 3) |
| `face_masks` | bool | (S, F) |

Padding: zeros, faces −1 (clamped to 0 when used, `MS:175`), and mask prefixes.

**`nice/`** (26 files, 3.6 GB, mode 0, all clouds full at 2048 points):

| files | shapes |
|---|---|
| ABC17 × 6 | 500 ×5 + 383 |
| ABC2 × 6 | 500 ×5 + 76 |
| ABC55 × 5 | 500 each |
| ABC58 × 5 | 500 each |
| blobs × 4 | 500 each |

- Total 12,459 shapes; P = 2048, k = 64, V = F = 2048.
- Some ABC shapes are degenerate (2–8 faces; e.g. `ABC17_train_mode0_5` shape 13 is an open quad).
  Their GWN sign is meaningless.

**`simulated/`** (20 files, 2.2 GB, mode −1, 9,998 shapes):
- 3,588 are padded (n < 2048, minimum 12 points; about 25 have n < 64).

**Other files:**
- `pretraining_blob.npz`: 1 shape, 2048 points, 1024 faces.
- `ABC17_validation_3.npz`, `ABC0_validation_0.npz`: see §9.6.
- `FundamentalFormPredictor_pretrained.pkl`: see §3.

**Normals of length 2** in `pretraining_blob.npz` and both validation files. The cause is
`normals.col(i) += (…).normalized()` with `+=` instead of `=` (`S3.cpp:500`, also `S3.cpp:665` in a
noisy path). The network normalizes, but the training-time fitting uses raw normals (§17).

**RAM.** About 350 MB per 500-shape file uncompressed, 75 % of it `neighbors`; phase 4 holds
≈ 16 GB.

## 11. Data generation (`MD`, `S3.cpp`)

**Blobs** (`MD:364-500`):
- `F = |x|² − (1 + α·2·pnoise3(scale·x, base = seed % 256))²` on 64³ over `±1.2(1+α)`.
- Parameters: α ∈ (0.1, 0.7), scale ∈ (0.5, 2), octaves 1–4, persistence 0.5, lacunarity 2.
- PyVista contour (winding flipped), then pyfqmr simplification to 2048 faces (aggressiveness 10).
- Kept if watertight; then `center_and_scale`. `n_blobs = 3000`, seed 42, split 0.8; 2000 used.
- Only 256 distinct noise fields (`base = seed % 256`).
- The `noise` package's Perlin is Ken Perlin's improved noise with that library's permutation
  table. Exact regeneration needs it; our own blobs only need the same distribution.

**ABC** (`MD:730-859`):
- STL from chunk .7z files, deduplicated, `center_and_scale`.
- Shuffle and split 0.8 (seed 57).
- Skipped if more than 2048 faces or vertices, if the orientation test fails, if not watertight,
  or if the cloud has fewer than 10 points.

**Mode 0 sampling:**
- 50/50 between `sample_uniform_point_cloud(2048)` (area-weighted, face normals,
  `S3.cpp:130-160, 596-676`) and FPS from 4 × 2048 uniform samples (start index 0,
  `S2.cpp:96-139`).
- No noise, no flips.

**Mode −1** (`sample_even_raycasted_point_cloud`, `S3.cpp:206-293`):
- 4 cameras at Fibonacci directions with a random global rotation (`S2.cpp:51-73`), at distance
  √3, looking at the centroid.
- A square grid of side 2√3 with spacing `U(0.01, 0.1)` per camera.
- First hit (FCPW), normals flipped to face the ray.
- Subsampled to at most 2048 by FPS or random shuffle, 50/50.

**Reproducibility.** OpenMP loops share an `mt19937` in several samplers (`S3.cpp:491`,
`S2.cpp:782`): data races, not reproducible. Our generators use per-item seeded RNGs.

## 12. Analytic test SDFs (`S2.{h,cpp}`, `shape_3d.py:541-1088`)

**Shapes.** 16 classes:

| shape | parameters |
|---|---|
| circle | r |
| pie | r, t |
| arc | sc, ra, rb |
| pill | a, b, r (the paper's "segment") |
| vesica | r, d |
| box | b |
| cross | b |
| isosceles triangle | q |
| triangle | v0..2 |
| hexagon | r |
| ellipse | ab |
| heart | — |
| pentagon | r |
| moon | d, ra, rb |
| trapezoid | ra, rb, he |
| quad | v0..3 |

The paper (App. E) uses 14: no heart, no isosceles.

**Parameter ranges** (`sample_SDF`, `shape_3d.py:964-1088`):
- `r_min = .1`, `r_max = .9`, `z ∈ [.1, .9]`, `t ∈ [π/6, 5π/6]`, as in App. E.
- Triangle and quad perturbations are uniform in [0, 1]², normalized, so they point into the first
  quadrant only. The paper says "random vectors".
- The same key is reused for several parameters (circle/hexagon/pentagon single key; quad reuses
  `subkey2` for h and the perturbation).

**Extrusion and revolution** (`S2.h:147-209`):
- Revolution about y: `φ₂(‖(x,z)‖ − o, y)`. The gradient is unguarded at `‖(x,z)‖ = 0`.
- Extrusion: `min(max(w), 0) + ‖max(w, 0)‖` with `w = (φ₂(x,y), |z| − h)`. This is exact outside
  and in the interior only when `φ₂` is exact.
- Paper: `o, h ∈ [0.1, 0.3]`, random choice between extrude and revolve.
- **The 70 test instances cannot be reconstructed**: no seed list, and `sample_SDF` is called only
  by the demo. Presumably 14 shapes × 5 instances. We draw our own with recorded seeds.

**Sampling** (`S2.cpp:322-803`):
- Ray-cast surface samples. All crossings along random lines through the shape's bounding sphere
  (step `max(|φ|, 1e-4)`) give area-*approximately*-uniform samples, with the normal from the
  analytic gradient.
- `sample_farthest_point_cloud` uses 8n candidates (here 8×, unlike the meshes' 4×).
- `sample_uneven_point_cloud` **swaps extrusion and revolution** (`S2.cpp:727, 736`) and uses
  `bbox_extrusion` for both. It is broken for revolved shapes and racy (OpenMP with a shared RNG).
- Fewer hits than requested cycles samples; 0 hits divides by zero.

## 13. The 2D variant (Fig. 10)

- Not implemented in the reference. `shape_2d.*` has the analytic 2D SDFs (§12), sphere tracing in
  2D (`S2.cpp:177-236`; ε = 1e-4, at most 1000 steps, t ≤ 4), and the extrusion/revolution
  wrappers. Nothing fits circles to a 2D point cloud.
- A "points as circles" method would need either a 2D network (new training data and a new
  training run) or a classical fit: a per-point osculating circle from a local quadratic
  `Q(s) = a0 + a1 s + a2 s²`.
  - The curvature is `κ = 2a2/(1 + a1²)^{3/2}` and the center `Q*(0) − sign/|κ| n*`.
  - The circle's SDF is `s(‖x − c‖ − 1/|κ|)`, blended with the same Eq. 25–26.
- The SPEC lists this as a user decision.

## 14. Demo and shader (`DM`, `SH`)

**Inputs.**
- `data/*.obj|*.ply` and the 16 analytic shapes. Meshes are `center_and_scale`d and uniformly
  sampled to 2048 points; analytic shapes use 1024 samples, 50/50 revolve/extrude, param
  `U(0.1, 0.3)` (`DM:400-480`).
- Model: the repository pkl; coefficients in 50k chunks; then `fit_tori_from_forms`
  (`DM:853-893`).
- Point areas from libigl (`DM:822-839`).

**GPU layout.** `GL_RGBA32F` textures, `texelFetch`; the shader maps `idx → (idx % W, idx / W)`.
- Points: 2 texels each, `(p.xyz, area)` and `(n.xyz, 0)`.
- Tori: `(c.xyz, 0)` and `(θ, φ, R, r)`.
- Point BVH (C++ SAH, 12 bins, `bvh.cpp:171-250`): 3 texels per node, primitive indices as float
  (exact below 2²⁴), a fixed `int stack[64]`.

**Per-pixel SDF** (`blendedTorusDistanceGradient`, `SH:1361-1485`): the same blend.
- `λ = lambdaScale·64/σ`, `σ = ½ max dist`.
- **Brute force over all points by default** (`neighborhood_size = −1`), with area weights.
- In accelerated mode (radius query, then k-NN with k ≤ 128 by scan) there are **no area weights**.
- The demo never passes `neighborhood_radius` to the render pass (it stays 0), so the radius
  branch finds nothing and falls back to kNN (`DM:321, 713, 1184-1231`).

**Sphere tracing.**
- Bounding sphere radius `0.8·|diag| + isovalue`.
- Step `t += |φ − iso|/‖∇φ‖`, with no Lipschitz factor and no guard for `‖∇φ‖ = 0`.
- Primary rays: 40 steps, ε = 1e-2. Shadow rays: 80 steps, ε = 1e-4. `t ≤ 20`.
- Normal = normalized analytic `∇φ`.

**Shading.**
- Lambert with light `(0, 10, 2)`, floor 0.3, plus `0.3·IBL-GGX` (metallic 0.07, roughness 0.1,
  F0 0.2) from the blurred cubemap `demo/assets/ForestBlurred{0..5}.png`, plus `0.6·Fresnel^80`.
- Then `sqrt` (gamma about 2).
- Colors: contour purple `(68, 74, 145)/255`, mesh tan, tori by golden-ratio HSV
  (`network.py:651-667`).
- Ground plane at `bbox_min.z − 1` but with normal +y (an axis mix-up), with a hard shadow.

**Camera.** At `(10 cos 5.4, 2.4, 10 sin 5.4)`, looking at the origin, focal 4.5,
`p = (2·coord − res)/(res.y·zoom)`.

**Cut plane.** Axis-aligned, optionally rotated.
- Color map over `±2`: white→blue `(0.0627, 0.5176, 0.9490)` with `pow(2t, 7)` inside,
  white→grapefruit `(0.9922, 0.4314, 0.3373)` outside.
- Isobands `0.9(1 − 0.1·smoothstep(−0.01, 0.01, sin(−60φ)))`.

**Bugs.**
- The true-gradient readout uses `.x` for z (`SH:~2297`).
- The quantity list has 10 entries but the enum has 9.
- The tori render is O(n) per pixel.

## 15. Normal orientation (Sec. 5.3, Fig. 27)

Not in the reference. The paper text (p. 14):

**Procedure.**
1. For each `p_i`, compute a weighted average `ñ_i` of `∇φ` evaluated at its neighbors `p_j`, with
   weights `exp(−M_i(p_i − p_j))` and
   `M_i(r) = (κmin²⟨r, v̂min⟩² + κmax²⟨r, v̂max⟩² + ε²⟨r, n_i⟩²)^{1/2}`.
2. Flip `n_i` if `⟨n_i, ñ_i⟩ < 0`.
3. Refit the tori and repeat until convergence. Fig. 27 shows 10 iterations, 30 % → 0.3 % flipped.

**Unspecified:**
- the neighbor set (we propose the same k = 64),
- `ε`,
- whether `∇φ` is evaluated with or without `p_j`'s own torus,
- the convergence test,
- whether all points flip simultaneously (Jacobi) or sequentially.

**Problem.** `∇φ` at a data point is evaluated with `r̂ = 0/0` for that point (§7.4). It needs the
guard `r̂ = 0` when `q = p_i`; then `p_i`'s term contributes `w ∇g_i` only.

## 16. Paper vs code (inference and method)

| # | topic | paper | code |
|---|---|---|---|
| 1 | evaluation neighborhood | `R_eval = 2C/λ`, fallback 32-NN | 32-NN always (`infer.py:12, 76`) |
| 2 | `λ` | constant `10³/D`, `D` = mean 64-NN distance | per query `64/σ_q`; `σ_q` over the 32 used points |
| 3 | `σ_x` | `½ max` over all of P (Eq. 26) | over the neighbors used (32 at inference, 64 in training) |
| 4 | Eq. 25 vs Eq. 1 | Eq. 25 | Eq. 25 (shifted); identical to Eq. 1 mathematically |
| 5 | area weights | omitted | omitted (option `use_areas`, never set in Python); the shader uses them in brute force |
| 6 | eigenvectors | Eq. 20 | the `(v1 ± disc, v2)` form (equivalent, §6.1) |
| 7 | degenerate curvature | — | 1e-10 clamps (JAX and C++ differ, §6.1) |
| 8 | 2D circles (Fig. 10) | illustrated | absent |
| 9 | orientation (Fig. 27) | described | absent |
| 10 | Laplacian (App. D) | formula | placeholder `Δg` (TODO) |
| 11 | analytic shapes | 14 | 16 |
| 12 | evaluation grid | 64³ mean abs error (Fig. 17–19), 128³ for Fig. 7 | no evaluation scripts in the repo |
| 13 | baselines | SHC, SSPD with Barill areas (libigl) | in `TD.cpp`, requires `get_point_areas()` first; FWN precompute bug |
| 14 | training | §9.7 | §9.7 |

## 17. Bugs and hazards (with consequences for the port)

1. **NaN sources at inference:**
   - `q = p_i` (r̂), `v = 0` (on the axis), `‖d‖ = 0` (core circle);
   - `κ = 0` (C++ without clamps);
   - `|u_z| > 1` in `acos`.
   No filter; a NaN torus poisons every query that uses it.
2. **`get_average_neighbor_distance`** data race (`TD.h:134-146`).
3. **`fast_winding_number`** uses empty `areas` unless `get_point_areas()` was called
   (`TD.cpp:1260-1314`).
4. **Raw vs normalized normal.**
   - Features and the basis use `n/‖n‖`; the fitting formulas use the raw `n` (C++ and JAX).
   - For non-unit input (PLY via open3d, the norm-2 training files) the torus is wrong: the center
     offset `a00 n` and `n*` are scaled, and `n*` is not unit.
   - Our API requires unit normals: validate with a tolerance and return an error. Normalizing is
     an explicit opt-in.
5. **Outliers in `get_neighbors`** leave −1 columns (§4).
6. **`sample_isosurface`** (`TD.cpp:841-881`) evaluates gradients at *uninitialized* sample
   locations before filling them, reseeds `random_device` per iteration, and has an unbounded loop.
7. **Training NaN gradients** from padded rows zero whole steps (§9.1). The fix in the port is to
   compute features and coefficients only for valid rows and mask attention.
8. **Duplicates and ties** make the center ambiguous (§4).
9. **`generate_validation_datasets`** is broken (wrong arguments, `MD:1442-1458`).
10. **Shader:** NaN on the axis in brute force (no guard), `.x` for `.z`, floor axis.
11. **`sample_uneven_point_cloud`** swaps extrusion and revolution (§12).
12. **`point_cloud_bounding_box_2d`** races on shared min/max (`S2.cpp:79-91`). The 3D version uses
    per-thread copies and is fine.

## 18. Precision

**Network.**
- f32 throughout (weights, features, attention, LayerNorm).
- JAX on a GPU may use TF32 for f32 matmuls by default precision. Coefficients computed on the
  authors' 3090 can therefore differ from CPU JAX by about 1e-3 relative.
- Our golden values for the network come from JAX on CPU (§20) and from the black-box build
  (`runs_tmp/pat-ref/`), with tolerances, never bitwise.

**Features.** `p/σ − p₀/σ` on absolute coordinates. For a cloud at distance `L` from the origin
with spacing `h`, the relative feature error is `~ε_f32·L/h`. That is about 1e-3 for `L/h = 10⁴`.
The difference-first form `(p − p₀)/σ` is exact to ε. Our default is difference-first in f64,
then cast; the reference's order is available as a compatibility mode (SPEC D-list).

**Fitting.** f64 in C++. The coefficient inputs are f32, so tori carry about 1e-7 relative error
regardless.

**Torus SDF.** `‖d‖ − |r|` with `|r| ~ 1/|κmax|`, and `R`, `c` ~ `1/|κmin|`.
- For nearly flat points, `1/|κmin|` reaches 10³–10⁸ in units where the cloud spans 2. The
  absolute error is then `ε·(R + r)`: about 1e-13 to 1e-8 in f64, but 1e-4 to 10 in f32.
- So **f32 evaluation is unsafe in flat regions** (the authors' GLSL demo lives with it).
- A stable rewrite evaluates relative to `o = Q*(0,0)`, which is on the torus.
  - `x = (q − o) + (o − c)` gives `v − R` and `x·a` as `O(‖q − o‖)` quantities computed from small
    numbers, where possible.
  - That needs `v − R = (v² − R²)/(v + R)`, with `v² − R²` expanded around `o`.
  - This is our addition (SPEC D-list), proposed for the f32 GPU tracer only.

**Blend weights.** `e^{±64}` fits in f32 (max `e^{88.7}`). The sum of 32 weights at most
`32·e^{64}` is fine. All points with brute force: `n·e^{64}` overflows f32 for `n > e^{24.7}`
(about 5·10¹⁰), so it is fine. Underflow of far weights is harmless.

**Training.** f32 with `NAN_EPSILON = 1e-10` guards, which are below f32 resolution near 1 but
still prevent division by 0.

## 19. Performance facts (paper, README)

**Precompute** (kd-tree + network + fitting):
- About 40 s for 1M points on a 3090 + Threadripper 3990X.
- Under 10 s for 100k, 1–4 min for 5–10M, over 10 min for 20M+.
- M3 Max: about 1 min for 50k, 8 min for 1M (README).
- Figure 20: 29M points in 12.5 min, about 4 ms per query.

**Queries.** 7.4e-5 s per query on the M3 Max (sequential); 1e-4 to 1e-3 s in general.
Fast winding number is 3.9e-6 s.

**FLOPs.** About 218 MFLOP per neighborhood forward (8 layers × ≈27 MFLOP; FFN 62 %, QKV/O 31 %,
attention 8 %). That gives 2.2e14 FLOP for 1M points: about 5 s at 45 TFLOPS (RTX 4080 f32 peak),
or about 60–110 s at 2–3.6 TFLOPS (this machine's 32-core Zen 2, AVX2, peak about 3.6 TFLOPS f32).

## 20. Survey scripts and measurements (scratchpad, not in the repo)

**Weights.** A restricted unpickler maps the 7 allowed globals to stubs or numpy, then writes npz
and safetensors (`pkl/inspect_pkl.py`, `pkl/to_st.py`). Comparing the two pkls gave 136 of 136
arrays differing.

**JAX CPU forward oracle.**
- `jaxref/forward.py`: the architecture of §2 re-declared with Flax 0.12.2 / JAX 0.8.2 on CPU,
  loaded from the npz. It does not import the reference, and so needs no `pat_bindings`.
- kNN by scipy `cKDTree(k = 64)`.
- Coefficients for `tooth` (2048), `snail` (4096) and `holey-bunny` (12,917).
- JAX CPU throughput is about 1,050 neighborhoods/s on this machine, under an outside load of about
  27 cores (the reference build agent's golden runs).

**Eigenvector identity check.** `jaxref/eig.py`.

**Framework prototypes** (SPEC §8). The candle 0.10.2 CPU forward on 10k holey-bunny
neighborhoods:
- It reproduces the JAX coefficients to a max abs error of 2.2e-5.
- First version, with `broadcast_matmul`, ran at 42 neighborhoods/s.
- Further numbers are in SPEC §8.2.
