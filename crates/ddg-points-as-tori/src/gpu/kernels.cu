// Points as Tori: our NVRTC kernels (SPEC.md §5.4, §7.5; NOTES-GPU.md).
//
// The blend, the torus SDF and the fit follow the reference (points-as-tori, MIT, Copyright (c)
// 2026 Nicole Feng) through our CPU code (`src/blend.rs`, `src/sdf.rs`, `src/fit.rs`), which they
// mirror operation by operation.
//
// Determinism (the workspace GPU rule): one thread per output, in the CPU's operation order;
// gathers only, no atomics, no reductions across threads. Compiled with --fmad=false (and IEEE
// division and square root), so every + − × / √ is the CPU's: kNN, fitting and the per-torus SDF
// are bitwise equal to the CPU (level 2). The blend's `exp` is CUDA's, not the host libm's; the
// blend is level 1 (bitwise run to run) and within a measured tolerance of the CPU (NOTES-GPU.md).

typedef unsigned int u32;
typedef unsigned long long u64;

#define LEAF 0xFFFFFFFFu
#define STACK 64

// ---------------------------------------------------------------------------------------------
// kd-tree (the CPU tree, flattened): nodes in pre-order, the left child follows its parent.
// ---------------------------------------------------------------------------------------------

struct Tree {
    const double* lo;     // 3 per node
    const double* hi;     // 3 per node
    const u32* start;
    const u32* end;
    const u32* right;     // LEAF for a leaf
    const double* pts;    // 3 per point, tree order
    const u32* idx;       // tree order -> input index
};

__device__ __forceinline__ double box_d2(const Tree& t, u32 ni, const double* q) {
    double s = 0.0;
    for (int a = 0; a < 3; ++a) {
        double lo = t.lo[3 * ni + a], hi = t.hi[3 * ni + a];
        double d = q[a] < lo ? lo - q[a] : (q[a] > hi ? q[a] - hi : 0.0);
        s += d * d;
    }
    return s;
}

__device__ __forceinline__ double dist2(const double* p, const double* q) {
    double dx = p[0] - q[0];
    double dy = p[1] - q[1];
    double dz = p[2] - q[2];
    return dx * dx + dy * dy + dz * dz;
}

__device__ __forceinline__ bool less_(double da, u32 ia, double db, u32 ib) {
    return da < db || (da == db && ia < ib);
}

// The k nearest points of q by (d², index), sorted, into (od, oi); returns the count (min(k, n)).
// Exact: the same set and order as `KdTree::knn` (the result is unique; pruning is strict).
template <int KMAX>
__device__ int knn(const Tree& t, u32 n, const double* q, int k, double* od, u32* oi) {
    if (k > (int)n) k = (int)n;
    int cnt = 0;
    double worst = 1.0 / 0.0;
    u32 st_n[STACK];
    double st_d[STACK];
    int sp = 1;
    st_n[0] = 0;
    st_d[0] = box_d2(t, 0, q);
    while (sp > 0) {
        --sp;
        u32 ni = st_n[sp];
        double bd = st_d[sp];
        if (cnt == k && bd > worst) continue;
        if (t.right[ni] == LEAF) {
            for (u32 s = t.start[ni]; s < t.end[ni]; ++s) {
                double d2 = dist2(&t.pts[3 * s], q);
                u32 id = t.idx[s];
                if (cnt < k) {
                    int pos = cnt++;
                    while (pos > 0 && less_(d2, id, od[pos - 1], oi[pos - 1])) {
                        od[pos] = od[pos - 1];
                        oi[pos] = oi[pos - 1];
                        --pos;
                    }
                    od[pos] = d2;
                    oi[pos] = id;
                } else if (less_(d2, id, worst, oi[k - 1])) {
                    int pos = k - 1;
                    while (pos > 0 && less_(d2, id, od[pos - 1], oi[pos - 1])) {
                        od[pos] = od[pos - 1];
                        oi[pos] = oi[pos - 1];
                        --pos;
                    }
                    od[pos] = d2;
                    oi[pos] = id;
                } else {
                    continue;
                }
                if (cnt == k) worst = od[k - 1];
            }
        } else {
            u32 l = ni + 1, r = t.right[ni];
            double dl = box_d2(t, l, q), dr = box_d2(t, r, q);
            u32 a, b;
            double da, db;
            if (dl <= dr) { a = r; da = dr; b = l; db = dl; } else { a = l; da = dl; b = r; db = dr; }
            if (!(cnt == k && da > worst)) { st_n[sp] = a; st_d[sp] = da; ++sp; }
            if (!(cnt == k && db > worst)) { st_n[sp] = b; st_d[sp] = db; ++sp; }
        }
    }
    return cnt;
}

// Number of points with d² ≤ r² (the CPU's `KdTree::radius` test).
__device__ u32 radius_count(const Tree& t, const double* q, double r2) {
    u32 st[STACK];
    int sp = 1;
    st[0] = 0;
    u32 c = 0;
    while (sp > 0) {
        u32 ni = st[--sp];
        if (box_d2(t, ni, q) > r2) continue;
        if (t.right[ni] == LEAF) {
            for (u32 s = t.start[ni]; s < t.end[ni]; ++s)
                if (dist2(&t.pts[3 * s], q) <= r2) ++c;
        } else {
            st[sp] = t.right[ni];
            st[sp + 1] = ni + 1;
            sp += 2;
        }
    }
    return c;
}

// The points with d² ≤ r², sorted by (d², index), into (od, oi) (capacity `cap` = their count).
__device__ void radius_fill(const Tree& t, const double* q, double r2, double* od, u32* oi, u32 cap) {
    u32 st[STACK];
    int sp = 1;
    st[0] = 0;
    u32 c = 0;
    while (sp > 0) {
        u32 ni = st[--sp];
        if (box_d2(t, ni, q) > r2) continue;
        if (t.right[ni] == LEAF) {
            for (u32 s = t.start[ni]; s < t.end[ni]; ++s) {
                double d2 = dist2(&t.pts[3 * s], q);
                if (d2 <= r2 && c < cap) { od[c] = d2; oi[c] = t.idx[s]; ++c; }
            }
        } else {
            st[sp] = t.right[ni];
            st[sp + 1] = ni + 1;
            sp += 2;
        }
    }
    // heap sort by (d², index): O(m log m), in place
    u32 m = c;
    for (u32 start = m / 2; start-- > 0;) {
        u32 root = start;
        for (;;) {
            u32 ch = 2 * root + 1;
            if (ch >= m) break;
            if (ch + 1 < m && less_(od[ch], oi[ch], od[ch + 1], oi[ch + 1])) ++ch;
            if (!less_(od[root], oi[root], od[ch], oi[ch])) break;
            double td = od[root]; od[root] = od[ch]; od[ch] = td;
            u32 ti = oi[root]; oi[root] = oi[ch]; oi[ch] = ti;
            root = ch;
        }
    }
    for (u32 end = m; end-- > 1;) {
        double td = od[0]; od[0] = od[end]; od[end] = td;
        u32 ti = oi[0]; oi[0] = oi[end]; oi[end] = ti;
        u32 root = 0;
        for (;;) {
            u32 ch = 2 * root + 1;
            if (ch >= end) break;
            if (ch + 1 < end && less_(od[ch], oi[ch], od[ch + 1], oi[ch + 1])) ++ch;
            if (!less_(od[root], oi[root], od[ch], oi[ch])) break;
            double t2 = od[root]; od[root] = od[ch]; od[ch] = t2;
            u32 t3 = oi[root]; oi[root] = oi[ch]; oi[ch] = t3;
            root = ch;
        }
    }
}

#define TREE_ARGS const double *lo, const double *hi, const u32 *start, const u32 *end, const u32 *right, const double *pts, const u32 *tidx
#define TREE_MAKE Tree t = {lo, hi, start, end, right, pts, tidx}

// kNN batch: `knn_batch` of the CPU (k ≤ KMAX), rows of min(k, n).
template <int KMAX>
__device__ void knn_batch_impl(u64 nq, const double* qs, u32 n, int k, TREE_ARGS, u32* out_i, double* out_d) {
    u64 q = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= nq) return;
    TREE_MAKE;
    double od[KMAX];
    u32 oi[KMAX];
    int c = knn<KMAX>(t, n, &qs[3 * q], k, od, oi);
    for (int j = 0; j < c; ++j) { out_i[q * c + j] = oi[j]; out_d[q * c + j] = od[j]; }
}

// Fitting neighborhoods: row i = [i] ++ (the k − 1 nearest j ≠ i), from the k nearest (as
// `fitting_neighborhoods`). Query points are the data points (input order).
template <int KMAX>
__device__ void fitting_impl(u32 n, const double* points, int k, TREE_ARGS, u32* rows) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    TREE_MAKE;
    double od[KMAX];
    u32 oi[KMAX];
    int c = knn<KMAX>(t, n, &points[3 * i], k, od, oi);
    u32* row = &rows[i * (u64)c];
    row[0] = (u32)i;
    int w = 1;
    for (int j = 0; j < c; ++j)
        if (oi[j] != (u32)i && w < c) row[w++] = oi[j];
}

#define KNN_ENTRY(K)                                                                                       \
    extern "C" __global__ void knn_batch_##K(u64 nq, const double* qs, u32 n, int k, TREE_ARGS, u32* oi,  \
                                             double* od) {                                                 \
        knn_batch_impl<K>(nq, qs, n, k, lo, hi, start, end, right, pts, tidx, oi, od);                     \
    }                                                                                                      \
    extern "C" __global__ void fitting_##K(u32 n, const double* points, int k, TREE_ARGS, u32* rows) {     \
        fitting_impl<K>(n, points, k, lo, hi, start, end, right, pts, tidx, rows);                         \
    }

KNN_ENTRY(32)
KNN_ENTRY(64)
KNN_ENTRY(128)

// ---------------------------------------------------------------------------------------------
// Torus SDF (src/sdf.rs) in f64. [formula: pat-torus-sdf] [formula: pat-torus-sdf-gradient]
// [formula: pat-torus-sdf-laplacian]
// ---------------------------------------------------------------------------------------------

struct TorusArrays {
    const double* center;  // 3 per torus
    const double* axis;    // 3 per torus
    const double* major;
    const double* minor;
};

#define SING_NONE 0
#define SING_AXIS 1
#define SING_CORE 2

__device__ __forceinline__ double sign_pos(double x) { return x >= 0.0 ? 1.0 : -1.0; }

// (rv, v, v − R, d2)
__device__ __forceinline__ void parts(const TorusArrays& T, u32 i, const double* x, double* rv, double* v, double* d1, double* d2) {
    const double* c = &T.center[3 * i];
    const double* u = &T.axis[3 * i];
    rv[0] = x[0] - c[0];
    rv[1] = x[1] - c[1];
    rv[2] = x[2] - c[2];
    double cr0 = rv[1] * u[2] - rv[2] * u[1];
    double cr1 = rv[2] * u[0] - rv[0] * u[2];
    double cr2 = rv[0] * u[1] - rv[1] * u[0];
    *v = sqrt(cr0 * cr0 + cr1 * cr1 + cr2 * cr2);
    *d2 = rv[0] * u[0] + rv[1] * u[1] + rv[2] * u[2];
    *d1 = *v - T.major[i];
}

__device__ __forceinline__ double torus_sdf(const TorusArrays& T, u32 i, const double* x) {
    double rv[3], v, d1, d2;
    parts(T, i, x, rv, &v, &d1, &d2);
    double m = T.minor[i];
    return sign_pos(m) * (sqrt(d1 * d1 + d2 * d2) - fabs(m));
}

// g, ∇g, Δg (when lap), singular kind
__device__ __forceinline__ int torus_sdf_grad(const TorusArrays& T, u32 i, const double* x, bool lap, double* g, double* gr, double* lg) {
    double rv[3], v, d1, d2;
    parts(T, i, x, rv, &v, &d1, &d2);
    const double* u = &T.axis[3 * i];
    double m = T.minor[i];
    double nd = sqrt(d1 * d1 + d2 * d2);
    double s = sign_pos(m);
    *g = s * (nd - fabs(m));
    *lg = 0.0;
    if (nd == 0.0) {
        gr[0] = 0.0; gr[1] = 0.0; gr[2] = 0.0;
        return SING_CORE;
    }
    double w1 = d1 / nd, w2 = d2 / nd;
    if (v == 0.0) {
        gr[0] = s * w2 * u[0];
        gr[1] = s * w2 * u[1];
        gr[2] = s * w2 * u[2];
        return SING_AXIS;
    }
    gr[0] = s * (w1 * (rv[0] - d2 * u[0]) / v + w2 * u[0]);
    gr[1] = s * (w1 * (rv[1] - d2 * u[1]) / v + w2 * u[1]);
    gr[2] = s * (w1 * (rv[2] - d2 * u[2]) / v + w2 * u[2]);
    if (lap) *lg = s * (1.0 / nd + d1 / (nd * v));
    return SING_NONE;
}

// ---------------------------------------------------------------------------------------------
// The blend (src/blend.rs) in f64. [formula: pat-blend-shift] [formula: pat-blend-weight-range]
// [formula: pat-blend-convex] [formula: pat-blend-gradient] [formula: pat-blend-gradient-exact]
// [formula: pat-blend-laplacian]
// ---------------------------------------------------------------------------------------------

#define RULE_PERQUERY 0
#define RULE_CONSTANT 1
#define RULE_PAPER 2

#define WANT_DIST 0
#define WANT_GRAD_FROZEN 1
#define WANT_GRAD_EXACT 2
#define WANT_LAP 3

#define NB_LIST 0   // a sorted list of indices
#define NB_ALL 1    // 0..n in index order

struct BlendParams {
    int rule;
    double c;             // PerQuery C
    double lambda_const;  // Constant λ or PaperHeuristic λ
    double isovalue;
    int want;
};

struct BlendOut {
    double distance;
    double grad[3];
    double lap;
    u32 axis, core, at_data;
    int degenerate;
};

__device__ __forceinline__ u32 nb_at(int kind, const u32* list, u32 j) { return kind == NB_ALL ? j : list[j]; }

__device__ void blend(const BlendParams& P, const double* pos, const TorusArrays& T, const double* area, const double* x,
                      int kind, const u32* list, u32 m, BlendOut* o) {
    // shift: σ_x = ½ max ρ, min ρ, argmax of (ρ², index)
    double m2 = -1.0 / 0.0, min2 = 1.0 / 0.0;
    u32 arg = 0;
    for (u32 t = 0; t < m; ++t) {
        u32 i = nb_at(kind, list, t);
        double r0 = x[0] - pos[3 * i], r1 = x[1] - pos[3 * i + 1], r2 = x[2] - pos[3 * i + 2];
        double d2 = r0 * r0 + r1 * r1 + r2 * r2;
        if (d2 > m2 || (d2 == m2 && i > nb_at(kind, list, arg))) { m2 = d2; arg = t; }
        min2 = fmin(min2, d2);
    }
    double sigma = 0.5 * sqrt(m2), rho_min = sqrt(min2);
    double lambda, shift;
    o->degenerate = 0;
    if (P.rule == RULE_PERQUERY) {
        if (sigma == 0.0) { o->degenerate = 1; return; }
        lambda = P.c / sigma;
        shift = sigma;
    } else {
        lambda = P.lambda_const;
        shift = rho_min;
    }
    o->axis = 0; o->core = 0; o->at_data = 0;
    double wsum = 0.0, gsum = 0.0;
    if (P.want == WANT_DIST) {
        for (u32 t = 0; t < m; ++t) {
            u32 i = nb_at(kind, list, t);
            double r0 = x[0] - pos[3 * i], r1 = x[1] - pos[3 * i + 1], r2 = x[2] - pos[3 * i + 2];
            double rho = sqrt(r0 * r0 + r1 * r1 + r2 * r2);
            double a = area ? area[i] : 1.0;
            double w = a * exp(-lambda * (rho - shift));
            gsum += w * torus_sdf(T, i, x);
            wsum += w;
        }
        o->distance = gsum / wsum - P.isovalue;
        o->grad[0] = 0.0; o->grad[1] = 0.0; o->grad[2] = 0.0;
        o->lap = 0.0;
        return;
    }
    bool lap = P.want == WANT_LAP;
    double asum[3] = {0.0, 0.0, 0.0}, nsum[3] = {0.0, 0.0, 0.0};
    double numer1 = 0.0, numer2 = 0.0, esum = 0.0, ersum = 0.0;
    for (u32 t = 0; t < m; ++t) {
        u32 i = nb_at(kind, list, t);
        double rv[3] = {x[0] - pos[3 * i], x[1] - pos[3 * i + 1], x[2] - pos[3 * i + 2]};
        double rho = sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]);
        double rh[3];
        if (rho > 0.0) {
            rh[0] = rv[0] / rho; rh[1] = rv[1] / rho; rh[2] = rv[2] / rho;
        } else {
            o->at_data += 1;
            rh[0] = 0.0; rh[1] = 0.0; rh[2] = 0.0;
        }
        double a = area ? area[i] : 1.0;
        double w = a * exp(-lambda * (rho - shift));
        double g, gr[3], lg;
        int sing = torus_sdf_grad(T, i, x, lap, &g, gr, &lg);
        if (sing == SING_AXIS) o->axis += 1;
        if (sing == SING_CORE) o->core += 1;
        for (int c = 0; c < 3; ++c) {
            asum[c] += w * (gr[c] - lambda * rh[c] * g);
            nsum[c] += w * rh[c];
        }
        gsum += w * g;
        wsum += w;
        esum += w * rho * g;
        ersum += w * rho;
        if (lap) {
            double inv = rho > 0.0 ? 2.0 / rho : 0.0;
            double dt = rh[0] * gr[0] + rh[1] * gr[1] + rh[2] * gr[2];
            numer1 += w * (lg - 2.0 * lambda * dt + g * lambda * (lambda - inv));
            numer2 += w * (inv - lambda);
        }
    }
    double A[3] = {asum[0] / wsum, asum[1] / wsum, asum[2] / wsum};
    double B[3] = {lambda * nsum[0] / wsum, lambda * nsum[1] / wsum, lambda * nsum[2] / wsum};
    double phi = gsum / wsum;
    double grad[3] = {A[0] + phi * B[0], A[1] + phi * B[1], A[2] + phi * B[2]};
    if (P.want == WANT_GRAD_EXACT && P.rule == RULE_PERQUERY && sigma > 0.0) {
        double dphi = -(esum - phi * ersum) / wsum;
        u32 j = nb_at(kind, list, arg);
        double far[3] = {x[0] - pos[3 * j], x[1] - pos[3 * j + 1], x[2] - pos[3 * j + 2]};
        double lf = sqrt(far[0] * far[0] + far[1] * far[1] + far[2] * far[2]);
        if (lf > 0.0) {
            double k = -(P.c / (sigma * sigma)) * 0.5 / lf * dphi;
            for (int c = 0; c < 3; ++c) grad[c] += k * far[c];
        }
    }
    double laplacian = 0.0;
    if (lap) {
        double first = numer1 / wsum + (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
        double second = (grad[0] * B[0] + grad[1] * B[1] + grad[2] * B[2]) +
                         phi * (lambda * numer2 / wsum + (B[0] * B[0] + B[1] * B[1] + B[2] * B[2]));
        laplacian = first + second;
    }
    o->distance = phi - P.isovalue;
    o->grad[0] = grad[0]; o->grad[1] = grad[1]; o->grad[2] = grad[2];
    o->lap = laplacian;
}

// Per query outputs: dist, grad (3), lap, flags (4 u32: axis, core, at_data, fallback | degenerate
// << 31), terms.
__device__ __forceinline__ void store(u64 q, const BlendOut& o, u32 fell_back, u32 terms, double* dist, double* grad, double* lapo,
                                      u32* flags) {
    dist[q] = o.distance;
    if (grad) { grad[3 * q] = o.grad[0]; grad[3 * q + 1] = o.grad[1]; grad[3 * q + 2] = o.grad[2]; }
    if (lapo) lapo[q] = o.lap;
    flags[5 * q] = o.axis;
    flags[5 * q + 1] = o.core;
    flags[5 * q + 2] = o.at_data;
    flags[5 * q + 3] = fell_back | ((u32)o.degenerate << 31);
    flags[5 * q + 4] = terms;
}

#define EVAL_ARGS                                                                                            \
    u64 nq, const double *qs, u32 n, int rule, double c, double lambda_const, double isovalue, int want,     \
        const double *pos, const double *center, const double *axisv, const double *major, const double *minor, \
        const double *area, TREE_ARGS, double *dist, double *grad, double *lapo, u32 *flags

#define EVAL_SETUP                                                                     \
    u64 q = (u64)blockIdx.x * blockDim.x + threadIdx.x;                                \
    if (q >= nq) return;                                                               \
    TREE_MAKE;                                                                         \
    TorusArrays T = {center, axisv, major, minor};                                     \
    BlendParams P = {rule, c, lambda_const, isovalue, want};                           \
    const double* x = &qs[3 * q];                                                      \
    BlendOut o;

// kNN neighborhoods (k < n, k ≤ KMAX).
template <int KMAX>
__device__ void eval_knn_impl(int k, EVAL_ARGS) {
    EVAL_SETUP;
    double od[KMAX];
    u32 oi[KMAX];
    int cnt = knn<KMAX>(t, n, x, k, od, oi);
    blend(P, pos, T, area, x, NB_LIST, oi, (u32)cnt, &o);
    store(q, o, 0, (u32)cnt, dist, grad, lapo, flags);
}

#define EVAL_ENTRY(K)                                                                                    \
    extern "C" __global__ void eval_knn_##K(int k, EVAL_ARGS) {                                          \
        eval_knn_impl<K>(k, nq, qs, n, rule, c, lambda_const, isovalue, want, pos, center, axisv, major,  \
                         minor, area, lo, hi, start, end, right, pts, tidx, dist, grad, lapo, flags);    \
    }

EVAL_ENTRY(32)
EVAL_ENTRY(64)
EVAL_ENTRY(128)

// All points in index order (`Neighborhood::All`, kNN with k ≥ n).
extern "C" __global__ void eval_all(EVAL_ARGS) {
    EVAL_SETUP;
    (void)t;
    blend(P, pos, T, area, x, NB_ALL, 0, n, &o);
    store(q, o, 0, n, dist, grad, lapo, flags);
}

// Radius neighborhoods, pass 1: the count of every query.
extern "C" __global__ void radius_counts(u64 nq, const double* qs, double r2, TREE_ARGS, u32* counts) {
    u64 q = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= nq) return;
    TREE_MAKE;
    counts[q] = radius_count(t, &qs[3 * q], r2);
}

// Radius neighborhoods, pass 2: fill each query's segment `[offsets[q], offsets[q] + counts[q])`
// of the scratch (sorted by (d², index)), blend over it; an empty set falls back to the
// `fallback_k` nearest (≤ 128).
extern "C" __global__ void eval_radius(double r2, int fallback_k, const u64* offsets, const u32* counts, double* sd, u32* si,
                                       EVAL_ARGS) {
    EVAL_SETUP;
    u32 m = counts[q];
    if (m > 0) {
        double* od = &sd[offsets[q]];
        u32* oi = &si[offsets[q]];
        radius_fill(t, x, r2, od, oi, m);
        blend(P, pos, T, area, x, NB_LIST, oi, m, &o);
        store(q, o, 0, m, dist, grad, lapo, flags);
    } else {
        double od[128];
        u32 oi[128];
        int cnt = knn<128>(t, n, x, fallback_k, od, oi);
        blend(P, pos, T, area, x, NB_LIST, oi, (u32)cnt, &o);
        store(q, o, 1, (u32)cnt, dist, grad, lapo, flags);
    }
}

// ---------------------------------------------------------------------------------------------
// F32Stable (SPEC.md §7.3, ours): the blend in f32 with the stable torus SDF. Coordinates are
// relative to the bounding-box center (`x32`, `pos32`); each torus is expanded around its data
// point p: with y = x − p, w = p − c (exact in f64 on the host, then rounded), z = y + w:
//   v² − v_p² = |y|² + 2 y·w − (y·u)² − 2 (y·u) d2p          (no O(R²) cancellation)
//   δ1 = (v² − v_p²)/(v + v_p),  δ2 = y·u,  d = d_p + δ
//   ‖d‖² − r² = K₂ + 2 d_p·δ + |δ|²,  g = s (‖d‖² − r²)/(‖d‖ + |r|)  (no O(r) cancellation)
// with v_p, d_p = (v_p − R, w·u), K₂ = ‖d_p‖² − r² from the host in f64. kNN stays in f64 (the
// neighbor sets are the CPU's). [formula: pat-torus-sdf-stable]
// ---------------------------------------------------------------------------------------------

struct Stable {
    const float* pos;   // 3 per point, relative to the origin
    const float* w;     // 3 per torus: p − c
    const float* axis;  // 3 per torus
    const float* vp;
    const float* d1p;
    const float* d2p;
    const float* k2;
    const float* minor;
    const float* area;
};

__device__ __forceinline__ int stable_sdf(const Stable& S, u32 i, const float* y, bool grad, bool lap, float* g, float* gr, float* lg) {
    const float* w = &S.w[3 * i];
    const float* u = &S.axis[3 * i];
    float yu = y[0] * u[0] + y[1] * u[1] + y[2] * u[2];
    float d2p = S.d2p[i], d1p = S.d1p[i], vp = S.vp[i];
    float A = (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) + 2.0f * (y[0] * w[0] + y[1] * w[1] + y[2] * w[2]) - yu * yu - 2.0f * yu * d2p;
    float zu = yu + d2p;
    float pz[3] = {y[0] + w[0] - zu * u[0], y[1] + w[1] - zu * u[1], y[2] + w[2] - zu * u[2]};
    float v = sqrtf(pz[0] * pz[0] + pz[1] * pz[1] + pz[2] * pz[2]);
    float den = v + vp;
    float e1 = den > 0.0f ? A / den : 0.0f;
    float e2 = yu;
    float d1 = d1p + e1, d2 = d2p + e2;
    float nd = sqrtf(d1 * d1 + d2 * d2);
    float m = S.minor[i];
    float s = m >= 0.0f ? 1.0f : -1.0f;
    float r = fabsf(m);
    float E = S.k2[i] + 2.0f * (d1p * e1 + d2p * e2) + (e1 * e1 + e2 * e2);
    float q = nd + r;
    *g = q > 0.0f ? s * (E / q) : 0.0f;
    *lg = 0.0f;
    if (!grad) return SING_NONE;
    if (nd == 0.0f) {
        gr[0] = 0.0f; gr[1] = 0.0f; gr[2] = 0.0f;
        return SING_CORE;
    }
    float w1 = d1 / nd, w2 = d2 / nd;
    if (v == 0.0f) {
        gr[0] = s * w2 * u[0]; gr[1] = s * w2 * u[1]; gr[2] = s * w2 * u[2];
        return SING_AXIS;
    }
    for (int c = 0; c < 3; ++c) gr[c] = s * (w1 * pz[c] / v + w2 * u[c]);
    if (lap) *lg = s * (1.0f / nd + d1 / (nd * v));
    return SING_NONE;
}

// The f32 blend in the cancellation-free arrangement (equal to App. D and to `blend` in exact
// arithmetic): pass 1 the shift, W and φ; pass 2 ∇φ = Σ w ∇g/W + λ Σ w r̂ (φ − g)/W (the App. D
// form A + φB cancels O(λ|g|) terms, ~1e3 here); pass 3 (Laplacian)
// Δφ = Σ w [Δg − 2λ r̂·(∇g − ∇φ) + (g − φ)(λ² − 2λ/ρ)]/W. The torus terms are recomputed per pass.
__device__ void blend32(const BlendParams& P, const Stable& S, const float* x, int kind, const u32* list, u32 m, BlendOut* o) {
    float m2 = -1.0f / 0.0f, min2 = 1.0f / 0.0f;
    u32 arg = 0;
    for (u32 t = 0; t < m; ++t) {
        u32 i = nb_at(kind, list, t);
        float r0 = x[0] - S.pos[3 * i], r1 = x[1] - S.pos[3 * i + 1], r2 = x[2] - S.pos[3 * i + 2];
        float d2 = r0 * r0 + r1 * r1 + r2 * r2;
        if (d2 > m2 || (d2 == m2 && i > nb_at(kind, list, arg))) { m2 = d2; arg = t; }
        min2 = fminf(min2, d2);
    }
    float sigma = 0.5f * sqrtf(m2), rho_min = sqrtf(min2);
    float lambda, shift;
    o->degenerate = 0;
    if (P.rule == RULE_PERQUERY) {
        if (sigma == 0.0f) { o->degenerate = 1; return; }
        lambda = (float)P.c / sigma;
        shift = sigma;
    } else {
        lambda = (float)P.lambda_const;
        shift = rho_min;
    }
    o->axis = 0; o->core = 0; o->at_data = 0;
    bool grad = P.want != WANT_DIST, lap = P.want == WANT_LAP;
    // pass 1: W, φ
    float wsum = 0.0f, gsum = 0.0f;
    for (u32 t = 0; t < m; ++t) {
        u32 i = nb_at(kind, list, t);
        float y[3] = {x[0] - S.pos[3 * i], x[1] - S.pos[3 * i + 1], x[2] - S.pos[3 * i + 2]};
        float rho = sqrtf(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
        float a = S.area ? S.area[i] : 1.0f;
        float w = a * expf(-lambda * (rho - shift));
        float g, gr[3], lg;
        stable_sdf(S, i, y, false, false, &g, gr, &lg);
        gsum += w * g;
        wsum += w;
    }
    float phi = gsum / wsum;
    o->distance = (double)phi - P.isovalue;
    o->grad[0] = 0.0; o->grad[1] = 0.0; o->grad[2] = 0.0;
    o->lap = 0.0;
    if (!grad) return;
    // pass 2: ∇φ
    float gw[3] = {0.0f, 0.0f, 0.0f}, nw[3] = {0.0f, 0.0f, 0.0f};
    float ew = 0.0f;
    for (u32 t = 0; t < m; ++t) {
        u32 i = nb_at(kind, list, t);
        float y[3] = {x[0] - S.pos[3 * i], x[1] - S.pos[3 * i + 1], x[2] - S.pos[3 * i + 2]};
        float rho = sqrtf(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
        float a = S.area ? S.area[i] : 1.0f;
        float w = a * expf(-lambda * (rho - shift));
        float g, gr[3], lg;
        int sing = stable_sdf(S, i, y, true, false, &g, gr, &lg);
        if (sing == SING_AXIS) o->axis += 1;
        if (sing == SING_CORE) o->core += 1;
        float rh[3] = {0.0f, 0.0f, 0.0f};
        if (rho > 0.0f) {
            rh[0] = y[0] / rho; rh[1] = y[1] / rho; rh[2] = y[2] / rho;
        } else {
            o->at_data += 1;
        }
        float dg = phi - g;
        for (int c = 0; c < 3; ++c) {
            gw[c] += w * gr[c];
            nw[c] += w * rh[c] * dg;
        }
        ew += w * rho * dg;
    }
    float gd[3];
    for (int c = 0; c < 3; ++c) gd[c] = gw[c] / wsum + lambda * (nw[c] / wsum);
    if (P.want == WANT_GRAD_EXACT && P.rule == RULE_PERQUERY && sigma > 0.0f) {
        // ∂φ/∂λ = −Σ w ρ (g − φ)/W = Σ w ρ (φ − g)/W
        float dphi = ew / wsum;
        u32 j = nb_at(kind, list, arg);
        float far[3] = {x[0] - S.pos[3 * j], x[1] - S.pos[3 * j + 1], x[2] - S.pos[3 * j + 2]};
        float lf = sqrtf(far[0] * far[0] + far[1] * far[1] + far[2] * far[2]);
        if (lf > 0.0f) {
            float k = -((float)P.c / (sigma * sigma)) * 0.5f / lf * dphi;
            for (int c = 0; c < 3; ++c) gd[c] += k * far[c];
        }
    }
    o->grad[0] = gd[0]; o->grad[1] = gd[1]; o->grad[2] = gd[2];
    if (!lap) return;
    // pass 3: Δφ (frozen λ)
    float lsum = 0.0f;
    for (u32 t = 0; t < m; ++t) {
        u32 i = nb_at(kind, list, t);
        float y[3] = {x[0] - S.pos[3 * i], x[1] - S.pos[3 * i + 1], x[2] - S.pos[3 * i + 2]};
        float rho = sqrtf(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
        float a = S.area ? S.area[i] : 1.0f;
        float w = a * expf(-lambda * (rho - shift));
        float g, gr[3], lg;
        stable_sdf(S, i, y, true, true, &g, gr, &lg);
        float rh[3] = {0.0f, 0.0f, 0.0f};
        float inv = 0.0f;
        if (rho > 0.0f) {
            rh[0] = y[0] / rho; rh[1] = y[1] / rho; rh[2] = y[2] / rho;
            inv = 2.0f / rho;
        }
        float dt = rh[0] * (gr[0] - gd[0]) + rh[1] * (gr[1] - gd[1]) + rh[2] * (gr[2] - gd[2]);
        lsum += w * (lg - 2.0f * lambda * dt + (g - phi) * lambda * (lambda - inv));
    }
    o->lap = (double)(lsum / wsum);
}

#define EVAL32_ARGS                                                                                       \
    u64 nq, const double *qs, u32 n, int rule, double c, double lambda_const, double isovalue, int want,  \
        const float *qs32, const float *pos32, const float *w32, const float *axis32, const float *vp,     \
        const float *d1p, const float *d2p, const float *k2, const float *minor32, const float *area32,    \
        TREE_ARGS, double *dist, double *grad, double *lapo, u32 *flags

#define EVAL32_SETUP                                                                   \
    u64 q = (u64)blockIdx.x * blockDim.x + threadIdx.x;                                \
    if (q >= nq) return;                                                               \
    TREE_MAKE;                                                                         \
    Stable S = {pos32, w32, axis32, vp, d1p, d2p, k2, minor32, area32};               \
    BlendParams P = {rule, c, lambda_const, isovalue, want};                           \
    const double* x = &qs[3 * q];                                                      \
    const float* x32 = &qs32[3 * q];                                                   \
    BlendOut o;

template <int KMAX>
__device__ void eval_knn_f32_impl(int k, EVAL32_ARGS) {
    EVAL32_SETUP;
    double od[KMAX];
    u32 oi[KMAX];
    int cnt = knn<KMAX>(t, n, x, k, od, oi);
    blend32(P, S, x32, NB_LIST, oi, (u32)cnt, &o);
    store(q, o, 0, (u32)cnt, dist, grad, lapo, flags);
}

#define EVAL32_ENTRY(K)                                                                                    \
    extern "C" __global__ void eval_knn_f32_##K(int k, EVAL32_ARGS) {                                      \
        eval_knn_f32_impl<K>(k, nq, qs, n, rule, c, lambda_const, isovalue, want, qs32, pos32, w32, axis32, \
                             vp, d1p, d2p, k2, minor32, area32, lo, hi, start, end, right, pts, tidx, dist,  \
                             grad, lapo, flags);                                                           \
    }

EVAL32_ENTRY(32)
EVAL32_ENTRY(64)
EVAL32_ENTRY(128)

extern "C" __global__ void eval_all_f32(EVAL32_ARGS) {
    EVAL32_SETUP;
    (void)t;
    (void)x;
    blend32(P, S, x32, NB_ALL, 0, n, &o);
    store(q, o, 0, n, dist, grad, lapo, flags);
}

extern "C" __global__ void eval_radius_f32(double r2, int fallback_k, const u64* offsets, const u32* counts, double* sd, u32* si,
                                           EVAL32_ARGS) {
    EVAL32_SETUP;
    u32 m = counts[q];
    if (m > 0) {
        double* od = &sd[offsets[q]];
        u32* oi = &si[offsets[q]];
        radius_fill(t, x, r2, od, oi, m);
        blend32(P, S, x32, NB_LIST, oi, m, &o);
        store(q, o, 0, m, dist, grad, lapo, flags);
    } else {
        double od[128];
        u32 oi[128];
        int cnt = knn<128>(t, n, x, fallback_k, od, oi);
        blend32(P, S, x32, NB_LIST, oi, (u32)cnt, &o);
        store(q, o, 1, (u32)cnt, dist, grad, lapo, flags);
    }
}

// ---------------------------------------------------------------------------------------------
// Torus fitting (src/fit.rs) in f64. [formula: pat-fundamental-forms] [formula: pat-mean-gauss]
// [formula: pat-principal-pair] [formula: pat-curvature-rayleigh] [formula: pat-principal-directions]
// [formula: pat-surface-normal] [formula: pat-torus-radii] [formula: pat-torus-sign]
// [formula: pat-torus-frame] [formula: pat-local-frame]
// ---------------------------------------------------------------------------------------------

__device__ __forceinline__ void cross3(const double* a, const double* b, double* o) {
    o[0] = a[1] * b[2] - a[2] * b[1];
    o[1] = a[2] * b[0] - a[0] * b[2];
    o[2] = a[0] * b[1] - a[1] * b[0];
}

__device__ __forceinline__ double sgn(double x) { return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0); }

__device__ __forceinline__ void torus_radii(double kmin, double kmax, double floor_, double* R, double* r) {
    double akmax = fmax(fabs(kmax), floor_);
    double akmin = fmax(fabs(kmin), floor_);
    *r = 1.0 / akmax;
    *R = 1.0 / akmin - sgn(kmax * kmin) * *r;
}

__device__ __forceinline__ double rayleigh(const double* ev, const double* I, const double* II, double eps, bool* cl) {
    double num = ev[0] * (II[0] * ev[0] + II[1] * ev[1]) + ev[1] * (II[1] * ev[0] + II[2] * ev[1]);
    double den = ev[0] * (I[0] * ev[0] + I[1] * ev[1]) + ev[1] * (I[1] * ev[0] + I[2] * ev[1]);
    *cl = fabs(den) < eps;
    den = sign_pos(den) * fmax(fabs(den), eps);
    return num / den;
}

__device__ __forceinline__ void pdir(const double* w, const double* s, const double* t, const double* n, double a10, double a01, double* o) {
    double ds[3] = {s[0] + a10 * n[0], s[1] + a10 * n[1], s[2] + a10 * n[2]};
    double dt[3] = {t[0] + a01 * n[0], t[1] + a01 * n[1], t[2] + a01 * n[2]};
    for (int c = 0; c < 3; ++c) o[c] = w[0] * ds[c] + w[1] * dt[c];
}

// One torus per point. coeffs: 6 f32 per point (a00, a01, a10, a11, a02, a20). Outputs: center
// (3), axis (3), major, minor (f64), clamp flags (bits: disc, eigvec, denominator, curvature,
// saddle, zero, eigen_switched).
extern "C" __global__ void fit_tori(u32 n, const double* P, const double* N, const float* A, double eps, double floor_, int geometric, int robust,
                                    double* center, double* axisv, double* major, double* minor, u32* clamps) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double* p = &P[3 * i];
    const double* nr = &N[3 * i];
    double a[6];
    for (int c = 0; c < 6; ++c) a[c] = (double)A[6 * i + c];
    double ln = sqrt(nr[0] * nr[0] + nr[1] * nr[1] + nr[2] * nr[2]);
    double cn[3] = {nr[0] / ln, nr[1] / ln, nr[2] / ln};
    // local frame of cn
    double uu[3] = {1.0, 0.0, 0.0};
    if (fabs(cn[0]) > 0.9) { uu[0] = 0.0; uu[1] = 1.0; }
    double s_hat[3], t_hat[3];
    cross3(cn, uu, s_hat);
    double ls = sqrt(s_hat[0] * s_hat[0] + s_hat[1] * s_hat[1] + s_hat[2] * s_hat[2]);
    s_hat[0] = s_hat[0] / ls; s_hat[1] = s_hat[1] / ls; s_hat[2] = s_hat[2] / ls;
    cross3(cn, s_hat, t_hat);
    double a00 = a[0], a01 = a[1], a10 = a[2], a11 = a[3], a02 = a[4], a20 = a[5];
    // fundamental forms
    double big_a = sqrt(1.0 + a01 * a01 + a10 * a10);
    double I[3] = {1.0 + a10 * a10, a10 * a01, 1.0 + a01 * a01};
    double II[3] = {2.0 * a20 / big_a, a11 / big_a, 2.0 * a02 / big_a};
    // principal pair terms
    double v1 = a20 * a01 * a01 - a02 * a10 * a10 - a02 + a20;
    double v2 = a11 * (a10 * a10 + 1.0) - 2.0 * a01 * a10 * a20;
    double term = a02 * (a10 * a10 + 1.0) - a01 * a10 * a11 + a20 * (a01 * a01 + 1.0);
    double disc2 = (1.0 + a01 * a01 + a10 * a10) * (a11 * a11 - 4.0 * a02 * a20) + term * term;
    u32 fl = 0;
    if (disc2 < eps) fl |= 1;
    double disc = sqrt(fmax(disc2, eps));
    double evp0 = v1 + disc;
    if (fabs(evp0) < eps) fl |= 2;
    evp0 = sign_pos(evp0) * fmax(fabs(evp0), eps);
    double evn1 = v2;
    if (fabs(evn1) < eps) fl |= 2;
    evn1 = sign_pos(evn1) * fmax(fabs(evn1), eps);
    double evp[2], evn[2];
    if (!robust) {
        evp[0] = evp0; evp[1] = v2;
        evn[0] = v1 - disc; evn[1] = evn1;
    } else {
        // EigenRule::Robust: the non-cancelling vector and its I-orthogonal complement J·I·e
        fl &= ~2u;
        double e0, e1, c0, c1, o0, o1, n0, n1;
        if (v1 >= 0.0) {
            e0 = v1 + disc; e1 = v2;
            c0 = -(I[1] * e0 + I[2] * e1); c1 = I[0] * e0 + I[1] * e1;
            evp[0] = e0; evp[1] = e1; evn[0] = c0; evn[1] = c1;
            o0 = v1 - disc; o1 = evn1; n0 = c0; n1 = c1;
        } else {
            e0 = v1 - disc; e1 = v2;
            c0 = -(I[1] * e0 + I[2] * e1); c1 = I[0] * e0 + I[1] * e1;
            evp[0] = c0; evp[1] = c1; evn[0] = e0; evn[1] = e1;
            o0 = evp0; o1 = v2; n0 = c0; n1 = c1;
        }
        double lo_ = sqrt(o0 * o0 + o1 * o1), ln_ = sqrt(n0 * n0 + n1 * n1);
        if (fabs((o0 * n1 - o1 * n0) / (lo_ * ln_)) > 1e-6) fl |= 64;
    }
    double lp = sqrt(evp[0] * evp[0] + evp[1] * evp[1]);
    double lq = sqrt(evn[0] * evn[0] + evn[1] * evn[1]);
    evp[0] = evp[0] / lp; evp[1] = evp[1] / lp;
    evn[0] = evn[0] / lq; evn[1] = evn[1] / lq;
    double vpos[3], vneg[3];
    pdir(evp, s_hat, t_hat, nr, a10, a01, vpos);
    pdir(evn, s_hat, t_hat, nr, a10, a01, vneg);
    double ns[3] = {(nr[0] - a10 * s_hat[0] - a01 * t_hat[0]) / big_a, (nr[1] - a10 * s_hat[1] - a01 * t_hat[1]) / big_a,
                    (nr[2] - a10 * s_hat[2] - a01 * t_hat[2]) / big_a};
    bool cp, cq;
    double kpos = rayleigh(evp, I, II, eps, &cp);
    double kneg = rayleigh(evn, I, II, eps, &cq);
    if (cp || cq) fl |= 4;
    double kmin, kmax;
    const double* vmin;
    if (fabs(kpos) < fabs(kneg)) { kmin = kpos; kmax = kneg; vmin = vpos; } else { kmin = kneg; kmax = kpos; vmin = vneg; }
    if (kmin == 0.0 || kmax == 0.0) fl |= 32;
    if (kmin * kmax < 0.0) fl |= 16;
    if (fabs(kmin) < floor_ || fabs(kmax) < floor_) fl |= 8;
    double R, r, s, side, akmin;
    if (!geometric) {
        torus_radii(kmin, kmax, floor_, &R, &r);
        s = (kmax < 0.0 && kmin < 0.0) || (kmax < 0.0 && kmin > 0.0) ? 1.0 : -1.0;
        side = s;
        akmin = fmax(fabs(kmin), floor_);
    } else {
        double km = fabs(kmin) >= floor_ ? kmin : (kmin < 0.0 ? -floor_ : floor_);
        double kx = fabs(kmax) >= floor_ ? kmax : (kmax < 0.0 ? -floor_ : floor_);
        torus_radii(km, kx, floor_, &R, &r);
        s = kx < 0.0 ? 1.0 : -1.0;
        side = km < 0.0 ? 1.0 : -1.0;
        akmin = fabs(km);
    }
    for (int c = 0; c < 3; ++c) center[3 * i + c] = p[c] + a00 * nr[c] - side / akmin * ns[c];
    double u[3];
    cross3(ns, vmin, u);
    double lu = fmax(sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]), eps);
    for (int c = 0; c < 3; ++c) axisv[3 * i + c] = u[c] / lu;
    major[i] = R;
    minor[i] = s * r;
    clamps[i] = fl;
}

// CUDA's exp, for the libm comparison of NOTES-GPU.md.
extern "C" __global__ void exp_values(u64 n, const double* x, double* y) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) y[i] = exp(x[i]);
}
