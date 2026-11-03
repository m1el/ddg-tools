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

// Near-field entries built on the device (level 1 only: pow() is CUDA's, not the CPU libm's, so the
// values may differ from the CPU's by an ulp). One CUDA block per near block q; the kernel of an
// element pair is evaluated for (min, max) as on the CPU (FaceKernels::near, VertexKernel::near), so
// mirrored blocks are exact transposes; the diagonal of a block (a, a) is 0 (never read).
__device__ __forceinline__ double dot3(const double* a, const double* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

extern "C" __global__ void near_fill_face(u32 n_blocks, const u32* __restrict__ block_row,
                                          const u32* __restrict__ near_col, const u64* __restrict__ near_off,
                                          const u32* __restrict__ start, const u32* __restrict__ end,
                                          const u32* __restrict__ order, const double* __restrict__ bary,
                                          const double* __restrict__ normal, double neg_s,
                                          double* __restrict__ near_val) {
    u32 q = blockIdx.x;
    if (q >= n_blocks) return;
    u32 a = block_row[q], b = near_col[q];
    u32 sa = start[a], ra = end[a] - sa, sb = start[b], cb = end[b] - sb;
    double* blk = near_val + 2 * near_off[q];
    for (u32 e = threadIdx.x; e < ra * cb; e += blockDim.x) {
        u32 r = e / cb, c = e % cb;
        double h0 = 0.0, h1 = 0.0;
        if (!(a == b && r == c)) {
            u32 i = order[sa + r], j = order[sb + c];
            u32 lo = i < j ? i : j, hi = i < j ? j : i;
            const double* xl = bary + 3 * (u64)lo;
            const double* xh = bary + 3 * (u64)hi;
            double d[3] = {xl[0] - xh[0], xl[1] - xh[1], xl[2] - xh[2]};
            double ns = dot3(normal + 3 * (u64)lo, d);
            double nt = dot3(normal + 3 * (u64)hi, d);
            double r2 = dot3(d, d);
            h0 = pow(r2, neg_s);
            h1 = 0.5 * (ns * ns + nt * nt) / (r2 * r2) * h0;
        }
        blk[2 * (u64)e] = h0;
        blk[2 * (u64)e + 1] = h1;
    }
}

extern "C" __global__ void near_fill_vertex(u32 n_blocks, const u32* __restrict__ block_row,
                                            const u32* __restrict__ near_col, const u64* __restrict__ near_off,
                                            const u32* __restrict__ start, const u32* __restrict__ end,
                                            const u32* __restrict__ order, const double* __restrict__ pos,
                                            double exponent, double* __restrict__ near_val) {
    u32 q = blockIdx.x;
    if (q >= n_blocks) return;
    u32 a = block_row[q], b = near_col[q];
    u32 sa = start[a], ra = end[a] - sa, sb = start[b], cb = end[b] - sb;
    double* blk = near_val + near_off[q];
    for (u32 e = threadIdx.x; e < ra * cb; e += blockDim.x) {
        u32 r = e / cb, c = e % cb;
        double h = 0.0;
        if (!(a == b && r == c)) {
            u32 i = order[sa + r], j = order[sb + c];
            u32 lo = i < j ? i : j, hi = i < j ? j : i;
            const double* xl = pos + 3 * (u64)lo;
            const double* xh = pos + 3 * (u64)hi;
            double d[3] = {xl[0] - xh[0], xl[1] - xh[1], xl[2] - xh[2]};
            h = pow(dot3(d, d), exponent);
        }
        blk[e] = h;
    }
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

// ---------------------------------------------------------------------------------------------
// Barnes-Hut energy and its exact gradient for frozen lists (src/energy/barnes_hut.rs,
// BhTree::face_partials / leaf_with_gradient / leaf_energy), split into gathers:
//   bh_own      one thread per face position: the own partials and the exact rows of its leaf;
//   bh_node     one thread per lumped entry: (b, u) of the entry;
//   bh_partner  one thread per partner slot (exact entry, target face): the partner partials;
//   bh_leaf_energy  one thread per leaf: its energy in the CPU's order;
//   bh_node_sums    one thread per node: (b_C, u_C) over the transposed lumped list;
//   bh_final    one thread per face position: own + partner (transposed exact list) + cluster part.
// The kernel K(X, N; Y) (src/kernel.rs) with the same operations (ipow by binary powering).

struct Exps {
    double alpha, beta;
    int alpha_int, beta_int;  // -1: not an integer (pow is used, not bitwise to the CPU)
};

__device__ __forceinline__ double ipow_(double x, u32 n) {
    double acc = 1.0;
    while (n > 0) {
        if (n & 1u) acc *= x;
        n >>= 1;
        if (n > 0) x *= x;
    }
    return acc;
}

__device__ __forceinline__ void nu_powers(const Exps& e, double nu, double* pa, double* sp1) {
    double a = fabs(nu);
    double p1 = e.alpha_int >= 0 ? ipow_(a, (u32)e.alpha_int - 1u) : pow(a, e.alpha - 1.0);
    *pa = p1 * a;
    *sp1 = copysign(p1, nu);
}

__device__ __forceinline__ double r_power(const Exps& e, double r2) {
    if (e.beta_int >= 0) {
        u32 k = (u32)e.beta_int;
        if (k % 2 == 0) return 1.0 / ipow_(r2, k / 2);
        return 1.0 / (ipow_(r2, k / 2) * sqrt(r2));
    }
    return pow(r2, -0.5 * e.beta);
}

__device__ __forceinline__ double k_value(const Exps& e, const double* v, const double* n) {
    double pa, sp1;
    nu_powers(e, dot3(n, v), &pa, &sp1);
    return pa * r_power(e, dot3(v, v));
}

// value, d_v, d_n of kernel_with_gradient_of_v.
__device__ __forceinline__ double k_grad(const Exps& e, const double* v, const double* n, double* dv, double* dn) {
    double r2 = dot3(v, v);
    double rb = r_power(e, r2);
    double pa, sp1;
    nu_powers(e, dot3(n, v), &pa, &sp1);
    double value = pa * rb;
    double c_n = e.alpha * sp1 * rb;
    double s = -e.beta * value / r2;
    for (int j = 0; j < 3; ++j) {
        dv[j] = n[j] * c_n + s * v[j];
        dn[j] = v[j] * c_n;
    }
    return value;
}

// Own partials (d/da, d/dX, d/dN) and the exact rows of the leaf of position ps.
extern "C" __global__ void bh_own(u32 np, Exps e, const u32* __restrict__ leaf_of, const u32* __restrict__ leaves,
                                  const u32* __restrict__ nstart, const u32* __restrict__ nend,
                                  const u32* __restrict__ lumped_off, const u32* __restrict__ lumped,
                                  const u32* __restrict__ exact_off, const u32* __restrict__ exact,
                                  const u32* __restrict__ row_off, const double* __restrict__ px,
                                  const double* __restrict__ pn, const double* __restrict__ pa,
                                  const double* __restrict__ cl_c, const double* __restrict__ cl_a,
                                  double* __restrict__ own, double* __restrict__ rows) {
    u32 ps = blockIdx.x * blockDim.x + threadIdx.x;
    if (ps >= np) return;
    u32 k = leaf_of[ps];
    u32 i = ps - nstart[leaves[k]];
    const double* xs = px + 3 * (u64)ps;
    const double* ns = pn + 3 * (u64)ps;
    double a_s = pa[ps];
    double o0 = 0.0, o1[3] = {0.0, 0.0, 0.0}, o2[3] = {0.0, 0.0, 0.0};
    double dv[3], dn[3], v[3];
    for (u32 q = lumped_off[k]; q < lumped_off[k + 1]; ++q) {
        u32 c = lumped[q];
        double ac = cl_a[c];
        const double* xc = cl_c + 3 * (u64)c;
        v[0] = xs[0] - xc[0];
        v[1] = xs[1] - xc[1];
        v[2] = xs[2] - xc[2];
        double val = k_grad(e, v, ns, dv, dn);
        o0 += val * ac;
        double s = a_s * ac;
        for (int j = 0; j < 3; ++j) {
            o1[j] = o1[j] + s * dv[j];
            o2[j] = o2[j] + s * dn[j];
        }
    }
    for (u32 q = exact_off[k]; q < exact_off[k + 1]; ++q) {
        u32 c = exact[q];
        double row = 0.0, dx[3] = {0.0, 0.0, 0.0}, dnn[3] = {0.0, 0.0, 0.0};
        for (u32 pt = nstart[c]; pt < nend[c]; ++pt) {
            if (pt == ps) continue;
            const double* xt = px + 3 * (u64)pt;
            double a_t = pa[pt];
            v[0] = xs[0] - xt[0];
            v[1] = xs[1] - xt[1];
            v[2] = xs[2] - xt[2];
            double val = k_grad(e, v, ns, dv, dn);
            row += a_t * val;
            for (int j = 0; j < 3; ++j) {
                dx[j] = dx[j] + a_t * dv[j];
                dnn[j] = dnn[j] + a_t * dn[j];
            }
        }
        rows[row_off[q] + i] = row;
        o0 += row;
        for (int j = 0; j < 3; ++j) {
            o1[j] = o1[j] + a_s * dx[j];
            o2[j] = o2[j] + a_s * dnn[j];
        }
    }
    double* o = own + 7 * (u64)ps;
    o[0] = o0;
    o[1] = o1[0];
    o[2] = o1[1];
    o[3] = o1[2];
    o[4] = o2[0];
    o[5] = o2[1];
    o[6] = o2[2];
}

// Value only (BhTree::leaf_energy): the exact rows of position ps.
extern "C" __global__ void bh_rows_value(u32 np, Exps e, const u32* __restrict__ leaf_of, const u32* __restrict__ leaves,
                                         const u32* __restrict__ nstart, const u32* __restrict__ nend,
                                         const u32* __restrict__ exact_off, const u32* __restrict__ exact,
                                         const u32* __restrict__ row_off, const double* __restrict__ px,
                                         const double* __restrict__ pn, const double* __restrict__ pa,
                                         double* __restrict__ rows) {
    u32 ps = blockIdx.x * blockDim.x + threadIdx.x;
    if (ps >= np) return;
    u32 k = leaf_of[ps];
    u32 i = ps - nstart[leaves[k]];
    const double* xs = px + 3 * (u64)ps;
    const double* ns = pn + 3 * (u64)ps;
    double v[3];
    for (u32 q = exact_off[k]; q < exact_off[k + 1]; ++q) {
        u32 c = exact[q];
        double row = 0.0;
        for (u32 pt = nstart[c]; pt < nend[c]; ++pt) {
            if (pt == ps) continue;
            const double* xt = px + 3 * (u64)pt;
            v[0] = xs[0] - xt[0];
            v[1] = xs[1] - xt[1];
            v[2] = xs[2] - xt[2];
            row += pa[pt] * k_value(e, v, ns);
        }
        rows[row_off[q] + i] = row;
    }
}

// (b, u) of every lumped entry; grad = 0: b only.
extern "C" __global__ void bh_node(u32 n_lumped, Exps e, u32 grad, const u32* __restrict__ lumped_leaf,
                                   const u32* __restrict__ lumped, const u32* __restrict__ leaves,
                                   const u32* __restrict__ nstart, const u32* __restrict__ nend,
                                   const double* __restrict__ px, const double* __restrict__ pn,
                                   const double* __restrict__ pa, const double* __restrict__ cl_c,
                                   double* __restrict__ node) {
    u32 q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= n_lumped) return;
    u32 l = leaves[lumped_leaf[q]];
    const double* xc = cl_c + 3 * (u64)lumped[q];
    double b = 0.0, u[3] = {0.0, 0.0, 0.0}, dv[3], dn[3], v[3];
    for (u32 ps = nstart[l]; ps < nend[l]; ++ps) {
        const double* xs = px + 3 * (u64)ps;
        double a_s = pa[ps];
        v[0] = xs[0] - xc[0];
        v[1] = xs[1] - xc[1];
        v[2] = xs[2] - xc[2];
        if (grad) {
            double val = k_grad(e, v, pn + 3 * (u64)ps, dv, dn);
            b += a_s * val;
            for (int j = 0; j < 3; ++j) u[j] = u[j] + -a_s * dv[j];
        } else {
            b += a_s * k_value(e, v, pn + 3 * (u64)ps);
        }
    }
    double* o = node + 4 * (u64)q;
    o[0] = b;
    o[1] = u[0];
    o[2] = u[1];
    o[3] = u[2];
}

// Partner partials (d/da_T, d/dX_T) of slot t = (exact entry q, target j).
extern "C" __global__ void bh_partner(u32 n_slots, Exps e, const u32* __restrict__ slot_entry,
                                      const u32* __restrict__ partner_off, const u32* __restrict__ exact_leaf,
                                      const u32* __restrict__ exact, const u32* __restrict__ leaves,
                                      const u32* __restrict__ nstart, const u32* __restrict__ nend,
                                      const double* __restrict__ px, const double* __restrict__ pn,
                                      const double* __restrict__ pa, double* __restrict__ partner) {
    u32 t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n_slots) return;
    u32 q = slot_entry[t];
    u32 pt = nstart[exact[q]] + (t - partner_off[q]);
    u32 l = leaves[exact_leaf[q]];
    const double* xt = px + 3 * (u64)pt;
    double a_t = pa[pt];
    double p0 = 0.0, p1[3] = {0.0, 0.0, 0.0}, dv[3], dn[3], v[3];
    for (u32 ps = nstart[l]; ps < nend[l]; ++ps) {
        if (ps == pt) continue;
        const double* xs = px + 3 * (u64)ps;
        double a_s = pa[ps];
        v[0] = xs[0] - xt[0];
        v[1] = xs[1] - xt[1];
        v[2] = xs[2] - xt[2];
        double val = k_grad(e, v, pn + 3 * (u64)ps, dv, dn);
        p0 += a_s * val;
        double s = -a_s * a_t;
        for (int j = 0; j < 3; ++j) p1[j] = p1[j] + s * dv[j];
    }
    double* o = partner + 4 * (u64)t;
    o[0] = p0;
    o[1] = p1[0];
    o[2] = p1[1];
    o[3] = p1[2];
}

// The energy of every leaf: lumped entries (b * A_C), then the exact rows (a_S * row).
extern "C" __global__ void bh_leaf_energy(u32 nl, const u32* __restrict__ leaves, const u32* __restrict__ nstart,
                                          const u32* __restrict__ nend, const u32* __restrict__ lumped_off,
                                          const u32* __restrict__ lumped, const u32* __restrict__ exact_off,
                                          const u32* __restrict__ row_off, const double* __restrict__ node,
                                          const double* __restrict__ cl_a, const double* __restrict__ pa,
                                          const double* __restrict__ rows, double* __restrict__ out) {
    u32 k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nl) return;
    u32 l = leaves[k];
    u32 s0 = nstart[l], ns = nend[l] - s0;
    double en = 0.0;
    for (u32 q = lumped_off[k]; q < lumped_off[k + 1]; ++q) en += node[4 * (u64)q] * cl_a[lumped[q]];
    for (u32 q = exact_off[k]; q < exact_off[k + 1]; ++q)
        for (u32 i = 0; i < ns; ++i) en += pa[s0 + i] * rows[row_off[q] + i];
    out[k] = en;
}

// (b_C, u_C) per node over the transposed lumped list (increasing entry).
extern "C" __global__ void bh_node_sums(u32 nn, const u32* __restrict__ lt_off, const u32* __restrict__ lt,
                                        const double* __restrict__ node, double* __restrict__ sums) {
    u32 n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= nn) return;
    double b = 0.0, u0 = 0.0, u1 = 0.0, u2 = 0.0;
    for (u32 i = lt_off[n]; i < lt_off[n + 1]; ++i) {
        const double* t = node + 4 * (u64)lt[i];
        b += t[0];
        u0 = u0 + 1.0 * t[1];
        u1 = u1 + 1.0 * t[2];
        u2 = u2 + 1.0 * t[3];
    }
    double* o = sums + 4 * (u64)n;
    o[0] = b;
    o[1] = u0;
    o[2] = u1;
    o[3] = u2;
}

// Per face: own + partner blocks (transposed exact list) + cluster part; out[7 f] by face f.
extern "C" __global__ void bh_final(u32 np, const u32* __restrict__ leaf_of, const u32* __restrict__ leaves,
                                    const u32* __restrict__ nstart, const u32* __restrict__ et_off,
                                    const u32* __restrict__ et, const u32* __restrict__ partner_off,
                                    const double* __restrict__ own, const double* __restrict__ partner,
                                    const double* __restrict__ beta, const double* __restrict__ ubig,
                                    const double* __restrict__ cl_c, const double* __restrict__ px,
                                    const double* __restrict__ pa, const u32* __restrict__ order,
                                    double* __restrict__ out) {
    u32 ps = blockIdx.x * blockDim.x + threadIdx.x;
    if (ps >= np) return;
    u32 k = leaf_of[ps];
    u32 c = leaves[k];
    u32 t = ps - nstart[c];
    const double* o = own + 7 * (u64)ps;
    double r0 = o[0], r1[3] = {o[1], o[2], o[3]};
    for (u32 i = et_off[k]; i < et_off[k + 1]; ++i) {
        const double* p = partner + 4 * ((u64)partner_off[et[i]] + t);
        r0 += p[0];
        r1[0] = r1[0] + 1.0 * p[1];
        r1[1] = r1[1] + 1.0 * p[2];
        r1[2] = r1[2] + 1.0 * p[3];
    }
    const double* xc = cl_c + 3 * (u64)c;
    const double* x = px + 3 * (u64)ps;
    const double* U = ubig + 3 * (u64)c;
    double d[3] = {x[0] - xc[0], x[1] - xc[1], x[2] - xc[2]};
    r0 += beta[c] + dot3(U, d);
    double a = pa[ps];
    double* w = out + 7 * (u64)order[ps];
    w[0] = r0;
    w[1] = r1[0] + a * U[0];
    w[2] = r1[1] + a * U[1];
    w[3] = r1[2] + a * U[2];
    w[4] = o[4];
    w[5] = o[5];
    w[6] = o[6];
}
