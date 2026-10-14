#!/usr/bin/env python3
"""Research diagnostics; no changes to production cone solvers.

Run with NumPy/SciPy, e.g. uv run --with numpy --with scipy python
scripts/investigate_camel.py --help. Build the Rust worker before UV experiments;
see crates/ddg-optimal-cones/experiments/camelshotnoise.md for commands and results.
All vertex numbers in printed diagnostics are zero-based.
"""
import argparse
import subprocess
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.optimize import minimize

DATA = Path(__file__).resolve().parents[1] / 'data/optimal-cones/optimal-cones-data/data'


def load(name, obj_path=None):
    p, uv, f, ft = [], [], [], []
    for line in (obj_path or DATA / f'{name}.obj').read_text().splitlines():
        w = line.split()
        if not w:
            continue
        if w[0] == 'v':
            p.append([float(x) for x in w[1:4]])
        elif w[0] == 'vt':
            uv.append([float(x) for x in w[1:3]])
        elif w[0] == 'f':
            assert len(w) == 4
            f.append([int(x.split('/')[0])-1 for x in w[1:]])
            ft.append([int(x.split('/')[1])-1 for x in w[1:]])
    cones = np.loadtxt(DATA / f'{name}-cones.txt')
    return np.array(p), np.array(f), np.array(uv)[np.array(ft)], cones[:, 0].astype(int)-1, cones[:, 1]


def operators(p, f, clamp=False):
    n = len(p)
    tri = p[f]
    cross = np.linalg.norm(np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]), axis=1)
    areas = cross / 2
    mass = np.bincount(f.ravel(), weights=np.repeat(areas/3, 3), minlength=n)
    mass /= areas.sum()
    angles = np.empty((len(f), 3))
    rows, cols, values = [], [], []
    for k in range(3):
        a, b = tri[:, (k+1)%3]-tri[:, k], tri[:, (k+2)%3]-tri[:, k]
        dot = np.einsum('ij,ij->i', a, b)
        angles[:, k] = np.arctan2(cross, dot)
        cot = dot / cross / 2
        i, j = f[:, (k+1)%3], f[:, (k+2)%3]
        rows.extend([i, j]); cols.extend([j, i]); values.extend([cot, cot])
    W = sparse.coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsc()
    if clamp:
        W.data = np.maximum(W.data, 0)
    L = sparse.diags(np.asarray(W.sum(axis=1)).ravel())-W
    edges = np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    edges, counts = np.unique(edges, axis=0, return_counts=True)
    bd = np.unique(edges[counts == 1])
    omega = 2*np.pi-np.bincount(f.ravel(), weights=angles.ravel(), minlength=n)
    omega[bd] = 0
    return L.tocsc(), mass, omega, bd, areas


def uv_measure(p, f, uv, areas):
    e1, e2 = p[f[:, 1]]-p[f[:, 0]], p[f[:, 2]]-p[f[:, 0]]
    a = np.linalg.norm(e1, axis=1)
    b = np.einsum('ij,ij->i', e1, e2)/a
    c = 2*areas/a
    U, V = uv[:, 1]-uv[:, 0], uv[:, 2]-uv[:, 0]
    J = np.stack([U/a[:, None], (V-U*(b/a)[:, None])/c[:, None]], axis=2)
    sv = np.linalg.svd(J, compute_uv=False)
    det = np.linalg.det(J)
    q = .5*np.log(np.abs(det))
    mean = np.average(q, weights=areas)
    return q, (q-mean)**2*areas/areas.sum(), sv[:, 0]/sv[:, 1], det < 0


def basis(p, f, sites, clamp=False):
    L, m, om, bd, areas = operators(p, f, clamp)
    ix = np.setdiff1d(np.arange(len(p)), bd)
    lu = splu(L[ix][:, ix])
    rhs = np.zeros((len(p), len(sites)))
    rhs[sites, np.arange(len(sites))] = 1
    G = np.zeros((len(p), len(sites)+len(bd)))
    G[ix, :len(sites)] = lu.solve(rhs[ix])
    G[ix, len(sites):] = lu.solve(-L[ix][:, bd].toarray())
    G[bd, len(sites):] = np.eye(len(bd))
    u0 = np.zeros(len(p)); u0[ix] = lu.solve(-om[ix])
    return u0, G, m, om, bd


def probe(name):
    p, f, uv, sites, angles = load(name)
    L, m, om, bd, areas = operators(p, f)
    tips = np.flatnonzero(om > 6)
    q, energy, anis, flips = uv_measure(p, f, uv, areas)
    spike_faces = np.any(np.isin(f, tips), axis=1)
    print(name, 'V,F,tips=', len(p), len(f), len(tips), 'area=', areas.sum())
    for label, mask in [('all', np.ones(len(f), dtype=bool)), ('spikes', spike_faces), ('rest', ~spike_faces), ('flipped', flips)]:
        if mask.any():
            print(label, 'faces', mask.sum(), 'area share', areas[mask].sum()/areas.sum(),
                  'Dstar contribution', energy[mask].sum(), 'flips', flips[mask].sum(),
                  'q quantiles', np.quantile(q[mask], [0, .5, 1]),
                  'anisotropy quantiles', np.quantile(anis[mask], [0, .5, .9, 1]))
    u0, G, _, _, _ = basis(p, f, sites)
    u = u0+G[:, :len(sites)]@angles
    print('linear E', .5*np.dot(m, u*u), 'u at tips', u[tips])
    print('tip diagonal', L.diagonal()[tips])
    if len(tips):
        W = L[tips].toarray()
        print('tip positive offdiag sum', np.maximum(W, 0).sum(axis=1)-L.diagonal()[tips])
    # A metric implied by exponentiating the Poisson state need not define valid triangles.
    lengths = np.stack([np.linalg.norm(p[f[:, k]]-p[f[:, (k+1)%3]], axis=1) for k in range(3)], axis=1)
    logs = np.log(lengths)+np.stack([.5*(u[f[:, k]]+u[f[:, (k+1)%3]]) for k in range(3)], axis=1)
    scaled = np.exp(logs-logs.max(axis=1)[:, None])
    invalid = 2*scaled.max(axis=1) >= scaled.sum(axis=1)
    print('invalid exponentiated-metric faces:', invalid.sum(), 'on spikes:', (invalid&spike_faces).sum())


def refit(name):
    p, f, uv, sites, angles = load(name)
    L, m, om, bd, areas = operators(p, f)
    tips = np.flatnonzero(om > 6)
    variants = [('original', p, False), ('positive-edge-cotan', p, True)]
    if len(tips):
        base = p.copy()
        for tip in tips:
            ring = np.unique(f[np.any(f == tip, axis=1)])
            ring = ring[ring != tip]
            base[tip] = p[ring].mean(axis=0)
        for height in [0., .0001, .001, .01, .1]:
            variants.append((f'needle-height-{height:g}', base+(p-base)*height, False))
    for variant, pos, clamp in variants:
        u0, G, m, _, _ = basis(pos, f, sites, clamp)
        root = np.sqrt(m)
        for free in [False, True]:
            H = G if free else G[:, :len(sites)]
            x = np.linalg.lstsq(root[:, None]*H, -root*u0, rcond=None)[0]
            print(variant, 'free' if free else 'fixed', 'max angle error', np.max(np.abs(x[:len(sites)]-angles)),
                  'E', .5*np.dot(m, (u0+H@x)**2), 'angles', x[:len(sites)], flush=True)
            if variant != 'original':
                continue
            losses = [('zero-tips', 0)] + [(loss, delta) for loss in ['huber', 'soft_l1', 'cauchy', 'asinh'] for delta in [.1, 1., 10.]]
            for loss, delta in losses:
                # Robustification is per log scale, with area outside rho(u).
                y = np.r_[angles, np.zeros(H.shape[1]-len(angles))]
                for it in range(100):
                    u = u0+H@y
                    if loss == 'zero-tips':
                        w = np.ones(len(p)); w[tips] = 0
                    elif loss == 'huber':
                        w = np.minimum(1, delta/np.maximum(np.abs(u), 1e-30))
                    elif loss == 'soft_l1':
                        w = 1/np.sqrt(1+(u/delta)**2)
                    elif loss == 'cauchy':
                        w = 1/(1+(u/delta)**2)
                    else:
                        z = u/delta
                        w = np.arcsinh(z)/(np.where(z == 0, 1., z)*np.sqrt(1+z*z))
                        w[z == 0] = 1
                    r = np.sqrt(m*w)
                    z = np.linalg.lstsq(r[:, None]*H, -r*u0, rcond=None)[0]
                    change = np.max(np.abs(z-y)); y = z
                    if change < 1e-8:
                        break
                print('  ', loss, delta, 'iters', it+1, 'change', change, 'max angle error', np.max(np.abs(y[:len(sites)]-angles)), 'angles', y[:len(sites)], flush=True)


class UVWorker:
    def __init__(self, name):
        binary = Path(__file__).resolve().parents[1] / 'target-camel-investigation/release/examples/camel_objective'
        self.proc = subprocess.Popen([str(binary), name], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)

    def request(self, command):
        line = self.raw(command)
        if line.startswith('error') or not line:
            return None
        return np.array([float(x) for x in line.split()])

    def raw(self, command):
        self.proc.stdin.write(command+'\n'); self.proc.stdin.flush()
        return self.proc.stdout.readline().strip()

    def measure(self, sites, angles):
        return self.request(' '.join(f'{v}:{a:.16e}' for v, a in zip(sites, angles)))

    def close(self):
        self.proc.stdin.close(); self.proc.wait()


def uv_opt(name):
    p, f, uv, sites, published = load(name)
    worker = UVWorker(name)
    print('columns: Dstar all absolute areas, Dstar excluding flipped, flips, spike flips, degenerate', flush=True)
    print('published UV', worker.request('published'), flush=True)
    print('BFF published angles', worker.measure(sites, published), flush=True)
    signs = np.sign(published)
    cap = np.abs(published).sum()
    for objective in [0, 1]:
        calls = 0
        def fun(x):
            nonlocal calls
            calls += 1
            metrics = worker.measure(sites, x)
            return 1e6 if metrics is None or metrics[4] else metrics[objective]
        # Same signs and total-angle budget; this is a fixed-position experiment.
        result = minimize(fun, published, method='SLSQP',
                          bounds=[(.01, 6.1) if s > 0 else (-6.1, -.01) for s in signs],
                          constraints=[{'type': 'ineq', 'fun': lambda x: cap-np.dot(signs, x)}],
                          options={'maxiter': 80, 'ftol': 1e-10, 'eps': 1e-4})
        print(json.dumps({'objective': ['all', 'exclude-flips'][objective], 'success': bool(result.success),
              'message': str(result.message), 'calls': calls, 'angles': result.x.tolist(),
              'max_angle_error': float(np.max(np.abs(result.x-published))),
              'sum_abs_angle': float(np.abs(result.x).sum()),
              'metrics': worker.measure(sites, result.x).tolist()}), flush=True)
    worker.close()


def proxy_scan(name):
    _, _, _, pub_sites, pub_angles = load(name)
    worker = UVWorker(name)
    print('published UV', worker.request('published'), flush=True)
    for boundary in ['fixed', 'free']:
        for lam in [.2, .1, .05, .02, .01, .005, .0045, .004, .0035, .003, .0025, .002, .001]:
            answer = worker.raw(f'generate {lam} {boundary}')
            if not answer.startswith('cones '):
                print(boundary, lam, answer, flush=True); continue
            pairs = [x.split(':') for x in answer.split()[1:]]
            sites = np.array([int(v) for v,a in pairs])
            angles = np.array([float(a) for v,a in pairs])
            metric = worker.measure(sites, angles)
            print(json.dumps({'boundary': boundary, 'lambda': lam, 'count': len(sites),
                  'sites': sites.tolist(), 'angles': angles.tolist(), 'total_angle': float(np.abs(angles).sum()),
                  'published_hits': len(set(sites)&set(pub_sites)),
                  'metrics': None if metric is None else metric.tolist()}), flush=True)
    worker.close()


def proxy_opt(name):
    _, _, _, pub_sites, published = load(name)
    worker = UVWorker(name)
    # Choose positions without consulting the published positions or angles.
    answer = worker.raw('generate 0.0045 free')
    assert answer.startswith('cones '), answer
    pairs = [x.split(':') for x in answer.split()[1:]]
    sites = np.array([int(v) for v,a in pairs])
    angles = np.array([float(a) for v,a in pairs])
    keep = np.argsort(-np.abs(angles))[:6]
    sites, angles = sites[keep], angles[keep]
    print('independent sites', sites, 'initial angles', angles, 'metrics', worker.measure(sites, angles), flush=True)
    cap = np.abs(published).sum()  # Budget is benchmark information, not locations.
    signs = np.sign(angles)
    calls = 0
    def fun(x):
        nonlocal calls
        calls += 1
        metrics = worker.measure(sites, x)
        if calls % 30 == 0:
            print('progress', calls, metrics, flush=True)
        return 1e6 if metrics is None or metrics[4] else metrics[0]
    result = minimize(fun, angles, method='SLSQP',
                      bounds=[(.01, 6.1) if s > 0 else (-6.1, -.01) for s in signs],
                      constraints=[{'type': 'ineq', 'fun': lambda x: cap-np.dot(signs, x)}],
                      options={'maxiter': 60, 'ftol': 1e-9, 'eps': 1e-4})
    print(json.dumps({'objective': 'all', 'success': bool(result.success), 'message': str(result.message),
          'calls': calls, 'sites': sites.tolist(), 'angles': result.x.tolist(),
          'published_hits': len(set(sites)&set(pub_sites)),
          'sum_abs_angle': float(np.abs(result.x).sum()),
          'metrics': worker.measure(sites, result.x).tolist()}), flush=True)
    worker.close()


def nonlinear(name):
    """Fixed-connectivity discrete-conformal Newton diagnostic, not a robust general solver."""
    p, f, uv, sites, angles = load(name)
    _, m, om, bd, areas = operators(p, f)
    ix = np.setdiff1d(np.arange(len(p)), bd)
    lengths = np.stack([np.linalg.norm(p[f[:, (k+1)%3]]-p[f[:, (k+2)%3]], axis=1) for k in range(3)], axis=1)
    target = np.zeros(len(p)); target[sites] = angles
    tips = np.flatnonzero(om > 6)

    def metric(s, hessian=True):
        logs = np.log(lengths)+np.stack([.5*(s[f[:, (k+1)%3]]+s[f[:, (k+2)%3]]) for k in range(3)], axis=1)
        top = logs.max(axis=1)
        ls = np.exp(logs-top[:, None])
        sorted_l = np.sort(ls, axis=1)
        c,b,a = sorted_l.T
        if np.any(c-(a-b) <= 0):
            return None
        four_area = np.sqrt((a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c)))
        cosine_num = np.stack([ls[:, (k+1)%3]**2+ls[:, (k+2)%3]**2-ls[:, k]**2 for k in range(3)], axis=1)
        theta = np.arctan2(four_area[:, None], cosine_num)
        K = 2*np.pi-np.bincount(f.ravel(), weights=theta.ravel(), minlength=len(p))
        log_area = np.log(four_area/4)+2*top
        if not hessian:
            return K, log_area
        rows, cols, vals = [], [], []
        for k in range(3):
            i,j = f[:, (k+1)%3],f[:, (k+2)%3]
            w = .5*cosine_num[:, k]/four_area
            rows.extend([i,j,i,j]); cols.extend([i,j,j,i]); vals.extend([w,w,-w,-w])
        L = sparse.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(len(p), len(p))).tocsc()
        return K, log_area, L

    s = np.zeros(len(p))
    direction = np.sin(np.arange(len(p)))
    direction[bd] = 0
    eps = 1e-7
    kplus = metric(eps*direction, False)[0]
    kminus = metric(-eps*direction, False)[0]
    jacobian = metric(s)[2]@direction
    print('initial curvature Jacobian directional-check relative error',
          np.linalg.norm(((kplus-kminus)/(2*eps)-jacobian)[ix])/np.linalg.norm(jacobian[ix]), flush=True)
    converged = False
    for it in range(100):
        K, log_area, L = metric(s)
        residual = K-target
        error = np.max(np.abs(residual[ix]))
        if it % 5 == 0 or error < 1e-8:
            q = .5*(log_area-np.log(areas)); mean = np.average(q, weights=areas)
            print('newton', it, 'max curvature residual', error, 'Dstar', np.dot(areas/areas.sum(), (q-mean)**2),
                  's range', (s.min(), s.max()), flush=True)
        if error < 1e-8:
            converged = True
            break
        delta = np.zeros(len(p))
        delta[ix] = splu(L[ix][:, ix]).solve(-residual[ix])
        alpha = min(1., 2./np.max(np.abs(delta)))
        norm = np.dot(residual[ix], residual[ix])
        for ls_iter in range(60):
            next_metric = metric(s+alpha*delta, False)
            if next_metric is not None:
                r = next_metric[0][ix]-target[ix]
                if np.dot(r,r) < (1.-1e-4*alpha)*norm:
                    break
            alpha *= .5
        else:
            print('line search failed', flush=True); break
        s += alpha*delta
        if alpha*np.max(np.abs(delta)) < 1e-10:
            print('stalled at triangle-inequality boundary; requires flips or a different model', flush=True); break
    print('CONVERGED' if converged else 'NOT CONVERGED', 'final maximum curvature residual',
          np.max(np.abs(metric(s, False)[0][ix]-target[ix])), flush=True)
    if len(tips):
        print('tip scales', s[tips], flush=True)


def needle(_name):
    """Exact regular six-triangle needle, flat disk after vertex scaling.

    Add unit area with zero scale to compare both objectives at fixed global scale.
    The boundary edge lengths of the fan are unchanged, so this is a valid
    discretely conformal metric with no apex curvature and positive triangle areas.
    It is not an angle-preserving affine map on each coarse triangle.
    """
    n, h = 6, .5
    print('h/rho, original area, Poisson u_tip, exact metric s_tip, linear D, area Dstar, flattened curvature')
    for ratio in [10., 100., 1000., 10000.]:
        rho = h/ratio
        slant = np.hypot(h, rho)
        base = 2*rho*np.sin(np.pi/n)
        altitude = np.sqrt(slant**2-base**2/4)
        area = n*base*altitude/2
        theta = 2*np.arcsin(base/(2*slant))
        omega = 2*np.pi-n*theta
        diagonal = n*np.tan(theta/2)
        u = -omega/diagonal
        s = 2*np.log(rho/slant)
        flattened_slant = np.exp(s/2)*slant
        flat_altitude = np.sqrt(flattened_slant**2-base**2/4)
        flat_area = n*base*flat_altitude/2
        q = .5*np.log(flat_area/area)
        Dlinear = (area/3)/(1+area)*u*u
        Darea = area/(1+area)**2*q*q
        curvature = 2*np.pi-2*n*np.arcsin(base/(2*flattened_slant))
        assert abs(curvature) < 1e-12
        print(ratio, area, u, s, Dlinear, Darea, curvature)


def export(name):
    root = Path(__file__).resolve().parents[1]
    log = root / 'runs_tmp/camel-objective-proxy-opt.log'
    record = [json.loads(line) for line in log.read_text().splitlines() if line.startswith('{')][-1]
    path = root / 'runs_tmp' / f'{name}-investigation.obj'
    worker = UVWorker(name)
    pairs = ' '.join(f'{v}:{a:.16e}' for v,a in zip(record['sites'], record['angles']))
    metrics = worker.request(f'save {path} {pairs}')
    worker.close()
    assert metrics is not None
    original_p, original_f, _, _, _ = load(name)
    p, f, uv, _, _ = load(name, path)
    assert np.array_equal(original_p, p), 'export changed original positions'
    assert np.array_equal(original_f, f), 'export changed original faces'
    areas = operators(p, f)[4]
    _, energy, _, flips = uv_measure(p, f, uv, areas)
    assert abs(energy.sum()-metrics[0]) < 1e-12
    assert flips.sum() == metrics[2]
    print('exported', path, 'verified unchanged geometry and independent area score', metrics)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['probe', 'refit', 'uv-opt', 'proxy-scan', 'proxy-opt', 'nonlinear', 'needle', 'export'])
    parser.add_argument('name', default='camelshotnoise', nargs='?')
    args = parser.parse_args()
    globals()[args.mode.replace('-', '_')](args.name)
