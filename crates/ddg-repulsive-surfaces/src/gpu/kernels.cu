// CUDA kernels of the GPU back end (NOTES-GPU.md), compiled at run time by NVRTC.
//
// Every kernel reproduces the operation order of the CPU code it replaces, one thread per output
// value, gathers only (no atomics, no shared-memory reductions): the results are bitwise identical
// from run to run whatever the device and launch configuration, and, compiled with --fmad=false
// (determinism level 2), bitwise equal to the CPU's (whose loops are compiled without FMA).
// The CPU counterparts are named in each comment (src/spatial/block_tree.rs, src/metric/hierarchical.rs).
//
// Layouts: element arrays are in BVH order ("position" p); a C-channel array is row-major [p * C + c];
// kernel coefficients are [entry * NK + k]; channel c uses kernel 0 if c < H, else kernel NK - 1
// (BlockTree::multiply_split).

typedef unsigned int u32;
typedef unsigned long long u64;

__device__ __forceinline__ u32 kidx(u32 c, u32 H, u32 NK) { return c < H ? 0u : NK - 1u; }

// BlockTree::multiply_core, step 1 (leaves): u = 0; u += w_p psi_p over the leaf's elements.
extern "C" __global__ void bct_up_leaves(u32 n_leaves, u32 C, const u32* __restrict__ leaves,
                                         const u32* __restrict__ start, const u32* __restrict__ end,
                                         const double* __restrict__ weight, const double* __restrict__ psi,
                                         double* __restrict__ up) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)n_leaves * C) return;
    u32 li = (u32)(t / C), c = (u32)(t % C);
    u32 a = leaves[li];
    double u = 0.0;
    for (u32 p = start[a]; p < end[a]; ++p) u += weight[p] * psi[(u64)p * C + c];
    up[(u64)a * C + c] = u;
}

// Step 1 (internal nodes of one depth, deepest first): u_a = u_left + u_right.
extern "C" __global__ void bct_up_internal(u32 count, u32 C, const u32* __restrict__ nodes,
                                           const u32* __restrict__ right, double* __restrict__ up) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)count * C) return;
    u32 a = nodes[t / C], c = (u32)(t % C);
    up[(u64)a * C + c] = up[(u64)(a + 1) * C + c] + up[(u64)right[a] * C + c];
}

// Step 2: the far field of every node, f = 0; f += k * u_col over its admissible blocks in order.
extern "C" __global__ void bct_far(u32 n_nodes, u32 C, u32 H, u32 NK, const u32* __restrict__ far_ptr,
                                   const u32* __restrict__ far_col, const double* __restrict__ far_coef,
                                   const double* __restrict__ up, double* __restrict__ far_own) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)n_nodes * C) return;
    u32 a = (u32)(t / C), c = (u32)(t % C);
    u32 k = kidx(c, H, NK);
    double f = 0.0;
    for (u32 q = far_ptr[a]; q < far_ptr[a + 1]; ++q)
        f += far_coef[(u64)q * NK + k] * up[(u64)far_col[q] * C + c];
    far_own[(u64)a * C + c] = f;
}

// Steps 3-4, one thread per (element row p, channel c): the downward accumulation as a gather along
// the chain root -> leaf (acc_child = own_child + acc_parent, the CPU's pre-order "acc[child] += acc[a]"),
// the far part w_p * g, then the near blocks of the chain in order, each row's columns in order,
// skipping the diagonal of a block (a, a).
extern "C" __global__ void bct_rows(u32 n, u32 C, u32 H, u32 NK, const u32* __restrict__ leaf_of,
                                    const u32* __restrict__ chain_ptr, const u32* __restrict__ chain,
                                    const u32* __restrict__ start, const u32* __restrict__ end,
                                    const u32* __restrict__ near_ptr, const u32* __restrict__ near_col,
                                    const u64* __restrict__ near_off, const double* __restrict__ near_val,
                                    const double* __restrict__ weight, const double* __restrict__ far_own,
                                    const double* __restrict__ psi, double* __restrict__ out) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)n * C) return;
    u32 p = (u32)(t / C), c = (u32)(t % C);
    u32 k = kidx(c, H, NK);
    u32 li = leaf_of[p];
    u32 c0 = chain_ptr[li], c1 = chain_ptr[li + 1];
    double g = far_own[(u64)chain[c0] * C + c];
    for (u32 i = c0 + 1; i < c1; ++i) g = far_own[(u64)chain[i] * C + c] + g;
    double acc = weight[p] * g;
    for (u32 i = c0; i < c1; ++i) {
        u32 a = chain[i];
        u32 row = p - start[a];
        for (u32 q = near_ptr[a]; q < near_ptr[a + 1]; ++q) {
            u32 b = near_col[q];
            u32 sb = start[b], cb = end[b] - start[b];
            const double* kr = near_val + (near_off[q] + (u64)row * cb) * NK + k;
            const double* x = psi + (u64)sb * C + c;
            if (a == b) {
                for (u32 j = 0; j < cb; ++j)
                    if (j != row) acc += kr[(u64)j * NK] * x[(u64)j * C];
            } else {
                for (u32 j = 0; j < cb; ++j) acc += kr[(u64)j * NK] * x[(u64)j * C];
            }
        }
    }
    out[t] = acc;
}

// FaceMetricTree::product, m = M x per face and psi = a m, in BVH order. Channels: 9 gradient
// channels per vector (9q + 3k + r) if hi, then 3 average channels per vector (off + 3q + k) if lo.
// xs: nb vectors of length 3|V|, one after the other.
extern "C" __global__ void face_m(u32 nf, u32 C, u32 off, u64 nd, const u32* __restrict__ order,
                                  const u32* __restrict__ faces, const double* __restrict__ df,
                                  const double* __restrict__ areas, const double* __restrict__ xs,
                                  double* __restrict__ m, double* __restrict__ psi) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)nf * C) return;
    u32 p = (u32)(t / C), ch = (u32)(t % C);
    u32 s = order[p];
    const u32* f = faces + 3 * (u64)s;
    double v;
    if (ch < off) {
        u32 q = ch / 9, rem = ch % 9, k = rem / 3, r = rem % 3;
        const double* x = xs + q * nd;
        const double* g = df + 9 * (u64)s;
        v = 0.0;
        for (u32 cc = 0; cc < 3; ++cc) v += x[3 * (u64)f[cc] + k] * g[3 * cc + r];
    } else {
        u32 j = ch - off, q = j / 3, k = j % 3;
        const double* x = xs + q * nd;
        v = (x[3 * (u64)f[0] + k] + x[3 * (u64)f[1] + k] + x[3 * (u64)f[2] + k]) / 3.0;
    }
    m[t] = v;
    psi[t] = v * areas[s];
}

// z = w (H~a m - H~ (a m)) + d m per face and channel, in BVH order (d = 0 without a diagonal).
extern "C" __global__ void face_z(u32 nf, u32 C, u32 off, double w_hi, double w_lo,
                                  const u32* __restrict__ order, const double* __restrict__ ha_perm,
                                  const double* __restrict__ d_hi, const double* __restrict__ d_lo,
                                  u32 has_d_hi, u32 has_d_lo, const double* __restrict__ m,
                                  double* __restrict__ hz) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)nf * C) return;
    u32 p = (u32)(t / C), ch = (u32)(t % C);
    u32 s = order[p];
    double w, d;
    u32 kern;
    if (ch < off) {
        w = w_hi;
        kern = 0;
        d = has_d_hi ? d_hi[s] : 0.0;
    } else {
        w = w_lo;
        kern = 1;
        d = has_d_lo ? d_lo[s] : 0.0;
    }
    double mm = m[t];
    hz[t] = w * (ha_perm[2 * (u64)p + kern] * mm - hz[t]) + d * mm;
}

// y_v = 2 sum_{S ni v} a_S (M^T row of S at v) . z_S, faces ascending (CSR incidence), per (v, q).
extern "C" __global__ void face_gather(u32 nv, u32 nb, u32 C, u32 off, u32 hi, u32 lo, u64 nd,
                                       const u32* __restrict__ inc_ptr, const u32* __restrict__ inc_face,
                                       const unsigned char* __restrict__ inc_corner,
                                       const double* __restrict__ areas, const double* __restrict__ df,
                                       const u32* __restrict__ position, const double* __restrict__ z,
                                       double* __restrict__ ys) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)nv * nb) return;
    u32 v = (u32)(t / nb), q = (u32)(t % nb);
    double acc0 = 0.0, acc1 = 0.0, acc2 = 0.0;
    for (u32 e = inc_ptr[v]; e < inc_ptr[v + 1]; ++e) {
        u32 s = inc_face[e];
        u32 c = inc_corner[e];
        double a = areas[s];
        const double* zs = z + (u64)position[s] * C;
        if (hi) {
            const double* g = df + 9 * (u64)s + 3 * c;
            const double* zh = zs + 9 * q;
            acc0 += a * (g[0] * zh[0] + g[1] * zh[1] + g[2] * zh[2]);
            acc1 += a * (g[0] * zh[3] + g[1] * zh[4] + g[2] * zh[5]);
            acc2 += a * (g[0] * zh[6] + g[1] * zh[7] + g[2] * zh[8]);
        }
        if (lo) {
            acc0 += a / 3.0 * zs[off + 3 * q];
            acc1 += a / 3.0 * zs[off + 3 * q + 1];
            acc2 += a / 3.0 * zs[off + 3 * q + 2];
        }
    }
    double* y = ys + q * nd + 3 * (u64)v;
    y[0] = 2.0 * acc0;
    y[1] = 2.0 * acc1;
    y[2] = 2.0 * acc2;
}

// VertexFractionalTree::product: psi[p][3q + k] = a_v x_q[3v + k], v = order[p].
extern "C" __global__ void vertex_psi(u32 n, u32 C, u64 nd, const u32* __restrict__ order,
                                      const double* __restrict__ mass, const double* __restrict__ xs,
                                      double* __restrict__ psi) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)n * C) return;
    u32 p = (u32)(t / C), ch = (u32)(t % C);
    u32 v = order[p], q = ch / 3, k = ch % 3;
    psi[t] = mass[v] * xs[q * nd + 3 * (u64)v + k];
}

// y_q[3v + k] = 2 a_v (H~a_v x_q[3v + k] - (H~ psi)[position v][3q + k]).
extern "C" __global__ void vertex_out(u32 n, u32 nb, u32 C, u64 nd, const u32* __restrict__ position,
                                      const double* __restrict__ mass, const double* __restrict__ ha,
                                      const double* __restrict__ xs, const double* __restrict__ hx,
                                      double* __restrict__ ys) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (u64)n * 3 * nb) return;
    u32 q = (u32)(t / (3 * (u64)n));
    u32 rem = (u32)(t % (3 * (u64)n));
    u32 v = rem / 3, k = rem % 3;
    double a = mass[v];
    ys[q * nd + rem] = 2.0 * a * (ha[v] * xs[q * nd + rem] - hx[(u64)position[v] * C + 3 * q + k]);
}

// The fixed-shape pairwise reduction (determinism, NOTES-GPU.md): out[i] = in[2i] + in[2i+1] over a
// zero-padded input of even length 2 * n_out (in[j] for j >= n_in is read as +0.0).
extern "C" __global__ void pairwise_level(u64 n_in, u64 n_out, const double* __restrict__ in,
                                          double* __restrict__ out) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_out) return;
    double a = 2 * i < n_in ? in[2 * i] : 0.0;
    double b = 2 * i + 1 < n_in ? in[2 * i + 1] : 0.0;
    out[i] = a + b;
}

// Element-wise products a_i b_i (the leaves of a dot product's pairwise tree).
extern "C" __global__ void mul_elems(u64 n, const double* __restrict__ a, const double* __restrict__ b,
                                     double* __restrict__ out) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = a[i] * b[i];
}
