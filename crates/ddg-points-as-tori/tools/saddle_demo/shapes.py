"""Saddle-dominated analytic shapes for the Points-as-Tori saddle demo (numpy only).

Conventions
-----------
* Principal curvatures are CONVEX POSITIVE with respect to the outward normal (a sphere of radius
  r has kappa = +1/r). The Monge coefficients handed to the torus fit use the fit's convention
  (height along the outward normal, so a convex patch has a20, a02 < 0).
* Monge patch of the fit (crate `fit.rs`, authors' `get_torus_parameters_python`):
  Q(s,t) = a00 + a10 s + a01 t + a11 s t + a20 s^2 + a02 t^2 along n, in the local frame
  (s_hat, t_hat) = local_frame(n); column order (a00, a01, a10, a11, a02, a20); physical units
  (already rescaled by sigma: a00 = sigma y0, slopes unscaled, curvature terms y/sigma).
  For an exact normal, a00 = a01 = a10 = 0 and the height is h = -1/2 x^T S x with S the shape
  operator (convex positive), so a20 = -S_ss/2, a02 = -S_tt/2, a11 = -S_st.
"""
import numpy as np

# ------------------------------------------------------------------ frame and coefficients


def local_frame(n):
    """(s_hat, t_hat) of the crate's `local_frame` and the authors' `orthonormal_basis` (rotation 0)."""
    n = np.asarray(n, float)
    u = np.where(np.abs(n[..., :1]) > 0.9, np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]))
    s = np.cross(n, u)
    s = s / np.linalg.norm(s, axis=-1, keepdims=True)
    t = np.cross(n, s)
    return s, t


def exact_coefficients(n, e1, k1, e2, k2):
    """Monge coefficients (a00, a01, a10, a11, a02, a20) of a surface point with unit outward normal
    n, orthonormal principal directions e1, e2 and convex-positive curvatures k1, k2."""
    s, t = local_frame(n)
    d = lambda a, b: np.sum(a * b, -1)
    S = lambda a, b: k1 * d(a, e1) * d(b, e1) + k2 * d(a, e2) * d(b, e2)
    z = np.zeros(len(n))
    return np.stack([z, z, z, -S(s, t), -0.5 * S(t, t), -0.5 * S(s, s)], -1)


def kappa_from_coefficients(c):
    """Principal curvatures of the Monge patch (generalized eigenvalues of (II, I)), returned
    CONVEX POSITIVE as (kmin, kmax) with |kmin| <= |kmax|."""
    a00, a01, a10, a11, a02, a20 = [c[:, i].astype(float) for i in range(6)]
    A = np.sqrt(1 + a01**2 + a10**2)
    I = np.stack([np.stack([1 + a10**2, a10 * a01], -1), np.stack([a10 * a01, 1 + a01**2], -1)], -2)
    II = np.stack([np.stack([2 * a20, a11], -1), np.stack([a11, 2 * a02], -1)], -2) / A[:, None, None]
    M = np.linalg.solve(I, II)
    ev = np.linalg.eigvals(M).real
    ev = -ev  # fit convention is convex negative
    order = np.argsort(np.abs(ev), axis=1)
    ev = np.take_along_axis(ev, order, 1)
    return ev[:, 0], ev[:, 1]


# ------------------------------------------------------------------ torus about z


class Torus:
    def __init__(self, R, r):
        self.R, self.r = R, r
        self.name = f'torus R {R} r {r}'

    def point(self, th, ph):
        R, r = self.R, self.r
        rho = R + r * np.cos(th)
        p = np.stack([rho * np.cos(ph), rho * np.sin(ph), r * np.sin(th)], -1)
        n = np.stack([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), np.sin(th)], -1)
        e_m = np.stack([-np.sin(th) * np.cos(ph), -np.sin(th) * np.sin(ph), np.cos(th)], -1)
        e_p = np.stack([-np.sin(ph), np.cos(ph), 0 * ph], -1)
        k_m = np.full_like(th, 1.0 / r)
        k_p = np.cos(th) / rho
        return p, n, e_m, k_m, e_p, k_p

    def sample(self, n, rng, cos_range=(-1.0, 1.0)):
        """Area-uniform samples with cos(theta) in cos_range (rejection on the area element)."""
        R, r = self.R, self.r
        th = np.empty(0)
        while len(th) < n:
            t = rng.uniform(0, 2 * np.pi, 4 * n)
            keep = (rng.uniform(0, 1, 4 * n) < (R + r * np.cos(t)) / (R + r)) & (np.cos(t) >= cos_range[0]) & (np.cos(t) <= cos_range[1])
            th = np.concatenate([th, t[keep]])
        th = th[:n]
        ph = rng.uniform(0, 2 * np.pi, n)
        return th, ph

    def sdf(self, x):
        rho = np.hypot(x[..., 0], x[..., 1])
        return np.hypot(rho - self.R, x[..., 2]) - self.r

    def nearest_param(self, x):
        rho = np.hypot(x[..., 0], x[..., 1])
        th = np.arctan2(x[..., 2], rho - self.R)
        ph = np.arctan2(x[..., 1], x[..., 0])
        return th, ph

    def gauss_at_nearest(self, x):
        th, _ = self.nearest_param(x)
        _, _, _, km, _, kp = self.point(th, np.zeros_like(th))
        return km * kp, np.cos(th)

    def saddle_area_fraction(self):
        return 0.5 - self.r / (np.pi * self.R)

    def mesh(self, nu=512, nv=256):
        th = np.linspace(0, 2 * np.pi, nv, endpoint=False)
        ph = np.linspace(0, 2 * np.pi, nu, endpoint=False)
        TH, PH = np.meshgrid(th, ph, indexing='ij')
        p = self.point(TH.ravel(), PH.ravel())[0]
        idx = lambda i, j: (i % nv) * nu + (j % nu)
        f = []
        for i in range(nv):
            for j in range(nu):
                a, b, c, d = idx(i, j), idx(i, j + 1), idx(i + 1, j + 1), idx(i + 1, j)
                f += [[a, b, c], [a, c, d]]
        return p, np.array(f, np.int64)


# ------------------------------------------------------------------ spool (hourglass) about z


class Spool:
    """Solid of revolution about z with profile, for u in [-pi/2, pi/2],
        z(u) = h sin u,   rho(u) = sqrt(a^2 + b^2 z^2) * cos u * sqrt(sum_{j<m} sin^{2j} u).
    The first factor is a one-sheet hyperboloid (rho^2 - b^2 z^2 = a^2: K < 0), the second equals
    sqrt(1 - sin^{2m} u) and closes the solid with round caps; rho is odd and z even about u = +-pi/2,
    so the surface is smooth (analytic) at the poles."""

    def __init__(self, h=0.9, a=0.22, b=0.55, m=8):
        self.h, self.a, self.b, self.m = h, a, b, m
        self.name = f'spool h {h} a {a} b {b} m {m}'
        self._table()

    def curve(self, u):
        s = np.sin(u)
        z = self.h * s
        cap = np.cos(u) * np.sqrt(sum(s ** (2 * j) for j in range(self.m)))
        rho = np.sqrt(self.a**2 + (self.b * z) ** 2) * cap
        return rho, z

    def derivs(self, u, eps=1e-4):
        """(rho, z), first and second derivatives by the 5-point stencil."""
        f = lambda v: np.stack(self.curve(v), -1)
        c0 = f(u)
        p1, m1, p2, m2 = f(u + eps), f(u - eps), f(u + 2 * eps), f(u - 2 * eps)
        d1 = (-p2 + 8 * p1 - 8 * m1 + m2) / (12 * eps)
        d2 = (-p2 + 16 * p1 - 30 * c0 + 16 * m1 - m2) / (12 * eps * eps)
        return c0, d1, d2

    def frame2d(self, u):
        """Profile point, outward unit normal (n_rho, n_z), meridian unit tangent and curvatures."""
        c, d1, d2 = self.derivs(u)
        sp = np.linalg.norm(d1, axis=-1)
        tan = d1 / sp[:, None]
        nrm = np.stack([d1[:, 1], -d1[:, 0]], -1) / sp[:, None]
        k_m = (d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]) / sp**3  # convex positive (CCW)
        k_p = nrm[:, 0] / c[:, 0]
        return c, nrm, tan, k_m, k_p

    def point(self, u, ph):
        c, nrm, tan, k_m, k_p = self.frame2d(u)
        cp, sp = np.cos(ph), np.sin(ph)
        p = np.stack([c[:, 0] * cp, c[:, 0] * sp, c[:, 1]], -1)
        n = np.stack([nrm[:, 0] * cp, nrm[:, 0] * sp, nrm[:, 1]], -1)
        e_m = np.stack([tan[:, 0] * cp, tan[:, 0] * sp, tan[:, 1]], -1)
        e_p = np.stack([-sp, cp, 0 * ph], -1)
        return p, n, e_m, k_m, e_p, k_p

    def _table(self, m=200001):
        u = np.linspace(-np.pi / 2, np.pi / 2, m)
        rho, z = self.curve(u)
        ds = np.hypot(np.diff(rho), np.diff(z))
        w = ds * 0.5 * (rho[1:] + rho[:-1])
        self.cdf = np.concatenate([[0], np.cumsum(w)])
        self.area = 2 * np.pi * self.cdf[-1]
        self.u_tab = u
        uu = np.linspace(-np.pi / 2, np.pi / 2, 8193)
        self.u_near = uu
        self.c_near = np.stack(self.curve(uu), -1)
        from scipy.spatial import cKDTree
        self.tree = cKDTree(self.c_near)

    def sample(self, n, rng, u_range=None):
        """Area-uniform samples (inverse CDF of rho ds on a fine table, then a uniform angle)."""
        lo, hi = 0.0, self.cdf[-1]
        if u_range is not None:
            lo, hi = np.interp(u_range, self.u_tab, self.cdf)
        q = rng.uniform(lo, hi, n)
        u = np.interp(q, self.cdf, self.u_tab)
        u = np.clip(u, -np.pi / 2 + 1e-6, np.pi / 2 - 1e-6)
        ph = rng.uniform(0, 2 * np.pi, n)
        return u, ph

    def nearest_u(self, rz):
        """Parameter of the nearest profile point to (rho, z) (table + Newton on (c - q).c' = 0)."""
        _, j = self.tree.query(rz)
        u = self.u_near[j].copy()
        for _ in range(12):
            c, d1, d2 = self.derivs(u, 1e-5)
            r = c - rz
            g = np.sum(r * d1, -1)
            hh = np.sum(d1 * d1, -1) + np.sum(r * d2, -1)
            step = np.where(np.abs(hh) > 1e-12, g / np.where(np.abs(hh) > 1e-12, hh, 1), 0)
            step = np.clip(step, -0.05, 0.05)
            u = np.clip(u - step, -np.pi / 2, np.pi / 2)
        return u

    def sdf(self, x, return_u=False):
        x = np.asarray(x, float).reshape(-1, 3)
        rz = np.stack([np.hypot(x[:, 0], x[:, 1]), x[:, 2]], -1)
        out = np.empty(len(x))
        us = np.empty(len(x))
        B = 200000
        for i in range(0, len(x), B):
            q = rz[i:i + B]
            u = self.nearest_u(q)
            c = np.stack(self.curve(u), -1)
            d = np.linalg.norm(c - q, axis=-1)
            inside = (np.abs(q[:, 1]) < self.h) & (q[:, 0] < self.curve(np.arcsin(np.clip(q[:, 1] / self.h, -1, 1)))[0])
            out[i:i + B] = np.where(inside, -d, d)
            us[i:i + B] = u
        return (out, us) if return_u else out

    def gauss_at_nearest(self, x):
        _, u = self.sdf(x, return_u=True)
        u = np.clip(u, -np.pi / 2 + 1e-6, np.pi / 2 - 1e-6)
        _, _, _, km, kp = self.frame2d(u)
        return km * kp, u

    def mesh(self, nu=512, nv=384):
        u = np.linspace(-np.pi / 2, np.pi / 2, nv + 1)[1:-1]
        ph = np.linspace(0, 2 * np.pi, nu, endpoint=False)
        rho, z = self.curve(u)
        V = [[0, 0, -self.h]]
        for i in range(len(u)):
            V += list(np.stack([rho[i] * np.cos(ph), rho[i] * np.sin(ph), np.full(nu, z[i])], -1))
        V.append([0, 0, self.h])
        V = np.array(V)
        top = len(V) - 1
        ring = lambda i, j: 1 + i * nu + (j % nu)
        F = []
        for j in range(nu):
            F.append([0, ring(0, j + 1), ring(0, j)])
            F.append([top, ring(len(u) - 1, j), ring(len(u) - 1, j + 1)])
        for i in range(len(u) - 1):
            for j in range(nu):
                a, b, c, d = ring(i, j), ring(i, j + 1), ring(i + 1, j + 1), ring(i + 1, j)
                F += [[a, b, c], [a, c, d]]
        F = np.array(F, np.int64)
        # orient outward (positive signed volume)
        vol = np.einsum('ij,ij->i', V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum() / 6
        if vol < 0:
            F = F[:, [0, 2, 1]]
        return V, F
