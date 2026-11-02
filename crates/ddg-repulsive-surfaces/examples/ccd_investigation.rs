//! Diagnostics of the CCD default (SPEC §7.3, §16 "CCD default (D20)"; NOTES-CCDINV.md).
//!
//! ```text
//! ccd_investigation run SCENE --steps N --threads T [--no-ccd] [--gap G] [--close-factor C]
//!                       [--csv PATH] [--obj PATH]
//!     Steps a scene (as `run_scene --no-stop`), checks every step end for intersections.
//!     With a locally instrumented session (env RS_CCDINV_DUMP=DIR, see NOTES-CCDINV.md) the
//!     step states are dumped; RS_CCDINV_CAP / RS_CCDINV_NO_GUARD vary the cap and the guard.
//! ccd_investigation check DIR --threads T [--l0 L]
//!     Exact swept test (cubic coplanarity roots, vertex-triangle and edge-edge) of every dumped
//!     step: start -> line-search trial, trial -> projection; plus the sheet distance and the
//!     largest remeshing displacement (the remesh pass is replayed from the dump).
//! ccd_investigation topo OBJ [OBJ ...]
//!     Genus, extents, area, volume; tree-cotree generators, the matrix V_ij = lk(g_i, g_j^+)
//!     (pushoff along the normal) and det(V - t V^T), an isotopy invariant of the embedding.
//! ccd_investigation klein SCENE --threads T
//!     CCD on an immersed input: touching / intersecting pairs at t = 0 and what the CCD reports.
//! ```

use std::collections::{HashMap, VecDeque};
use std::io::Write as _;
use std::path::PathBuf;

use ddg_repulsive_surfaces::ccd::{self, CcdOptions};
use ddg_repulsive_surfaces::remesh::{RemeshOptions, remesh};
use ddg_repulsive_surfaces::step::{CSV_HEADER, StopRule, csv_line, scene_options};
use ddg_repulsive_surfaces::{SceneSpec, Session, Surface};
use rayon::prelude::*;

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scl(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
fn lerp(a: V3, b: V3, t: f64) -> V3 {
    add(a, scl(sub(b, a), t))
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("run") => cmd_run(&args[1..]),
        Some("check") => cmd_check(&args[1..]),
        Some("topo") => cmd_topo(&args[1..]),
        Some("klein") => cmd_klein(&args[1..]),
        _ => {
            eprintln!("usage: ccd_investigation run|check|topo|klein ...");
            std::process::exit(2);
        }
    }
}

fn opt<T: std::str::FromStr>(args: &[String], name: &str) -> Option<T> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .map(|s| s.parse().ok().unwrap_or_else(|| panic!("bad value for {name}")))
}
fn flag(args: &[String], name: &str) -> bool {
    args.iter().any(|a| a == name)
}

// ---------------------------------------------------------------------------------------------
// run
// ---------------------------------------------------------------------------------------------

fn cmd_run(args: &[String]) {
    let scene = &args[0];
    let spec = SceneSpec::load(scene).expect("scene");
    let mut o = scene_options(&spec).expect("options");
    o.max_iterations = opt(args, "--steps").unwrap_or(200);
    o.threads = opt(args, "--threads");
    o.stop = StopRule::Off;
    if flag(args, "--no-ccd") {
        o.line_search.ccd = None;
    } else {
        let mut c = o.line_search.ccd.unwrap_or_default();
        if let Some(g) = opt(args, "--gap") {
            c.gap = g;
        }
        if let Some(cf) = opt(args, "--close-factor") {
            c.close_factor = cf;
        }
        o.line_search.ccd = Some(c);
    }
    let mut csv = opt::<PathBuf>(args, "--csv").map(|p| {
        let mut w = std::io::BufWriter::new(std::fs::File::create(p).expect("csv"));
        writeln!(w, "{CSV_HEADER},intersecting_pairs,ccd_contact").unwrap();
        w
    });
    let mut s = Session::new(spec.surface.clone(), spec.positions.clone(), o.clone()).expect("session");
    let e0 = s.energy();
    let mut total = 0.0;
    let mut n_int_steps = 0;
    let mut capped = 0;
    while !s.status().is_terminal() {
        let r = match s.step() {
            Ok(r) => r,
            Err(e) => {
                eprintln!("step error: {e}");
                break;
            }
        };
        total += r.seconds.total;
        let ni = ccd::intersecting_pairs(s.surface(), s.positions()).map_or(usize::MAX, |v| v.len());
        if ni > 0 {
            n_int_steps += 1;
        }
        let contact = r.ccd.as_ref().and_then(|c| c.contact).map_or(String::new(), |c| format!("{:?}", c.primitive));
        if r.ccd.as_ref().is_some_and(|c| c.contact.is_some()) {
            capped += 1;
        }
        if let Some(w) = csv.as_mut() {
            writeln!(w, "{},{},\"{}\"", csv_line(&r), ni, contact).unwrap();
            w.flush().unwrap();
        }
        if r.iteration % 20 == 0 {
            eprintln!(
                "step {} F={:.6e} tau={:.3e} tau0={:.3e} tau_ccd={:?} faces={} int={}",
                r.iteration, r.energy_after, r.tau, r.tau0, r.tau_ccd, r.faces, ni
            );
        }
    }
    let n = s.iteration().max(1) as f64;
    if let Some(p) = opt::<PathBuf>(args, "--obj") {
        ddg_repulsive_surfaces::write_obj(&p, s.surface(), s.positions()).expect("obj");
    }
    println!(
        "summary scene={scene} steps={} status={:?} energy0={e0:.6e} energy={:.6e} ms_per_step={:.1} \
         ccd_capped_steps={capped} steps_with_end_intersections={n_int_steps} faces={}",
        s.iteration(),
        s.status(),
        s.energy(),
        total * 1e3 / n,
        s.surface().n_faces()
    );
}

// ---------------------------------------------------------------------------------------------
// dumps
// ---------------------------------------------------------------------------------------------

struct Dump {
    iteration: u64,
    accepted: bool,
    tau: f64,
    faces0: Vec<[u32; 3]>,
    p0: Vec<V3>,
    moved: Vec<V3>,
    proj: Vec<V3>,
    faces1: Vec<[u32; 3]>,
    p1: Vec<V3>,
}

fn read_dump(path: &std::path::Path) -> Dump {
    let b = std::fs::read(path).expect("dump");
    let mut at = 0usize;
    let mut u = || {
        let x = u64::from_le_bytes(b[at..at + 8].try_into().unwrap());
        at += 8;
        x
    };
    let iteration = u();
    let accepted = u() == 1;
    let tau = f64::from_bits(u());
    let nv0 = u() as usize;
    let nf0 = u() as usize;
    let nv1 = u() as usize;
    let nf1 = u() as usize;
    let mut at = 56usize;
    let faces = |n: usize, at: &mut usize| -> Vec<[u32; 3]> {
        (0..n)
            .map(|_| {
                let f = [0, 1, 2].map(|k| u32::from_le_bytes(b[*at + 4 * k..*at + 4 * k + 4].try_into().unwrap()));
                *at += 12;
                f
            })
            .collect()
    };
    let faces0 = faces(nf0, &mut at);
    let pts = |n: usize, at: &mut usize| -> Vec<V3> {
        (0..n)
            .map(|_| {
                let p = [0, 1, 2].map(|k| f64::from_le_bytes(b[*at + 8 * k..*at + 8 * k + 8].try_into().unwrap()));
                *at += 24;
                p
            })
            .collect()
    };
    let p0 = pts(nv0, &mut at);
    let moved = pts(nv0, &mut at);
    let proj = pts(nv0, &mut at);
    let faces1 = {
        let f = |n: usize, at: &mut usize| -> Vec<[u32; 3]> {
            (0..n)
                .map(|_| {
                    let f = [0, 1, 2].map(|k| u32::from_le_bytes(b[*at + 4 * k..*at + 4 * k + 4].try_into().unwrap()));
                    *at += 12;
                    f
                })
                .collect()
        };
        f(nf1, &mut at)
    };
    let p1 = pts(nv1, &mut at);
    Dump {
        iteration,
        accepted,
        tau,
        faces0,
        p0,
        moved,
        proj,
        faces1,
        p1,
    }
}

// ---------------------------------------------------------------------------------------------
// exact swept test
// ---------------------------------------------------------------------------------------------

/// Roots in [0, 1] of the cubic through the values at 0, 1/3, 2/3, 1.
fn cubic_roots(f: impl Fn(f64) -> f64) -> Vec<f64> {
    let y = [f(0.0), f(1.0 / 3.0), f(2.0 / 3.0), f(1.0)];
    // monomial coefficients c0 + c1 t + c2 t^2 + c3 t^3 (Newton form on equispaced nodes)
    let h = 1.0 / 3.0;
    let d1 = [(y[1] - y[0]) / h, (y[2] - y[1]) / h, (y[3] - y[2]) / h];
    let d2 = [(d1[1] - d1[0]) / (2.0 * h), (d1[2] - d1[1]) / (2.0 * h)];
    let d3 = (d2[1] - d2[0]) / (3.0 * h);
    // p(t) = y0 + d1_0 t + d2_0 t (t - h) + d3 t (t - h)(t - 2h)
    let c3 = d3;
    let c2 = d2[0] - 3.0 * h * d3;
    let c1 = d1[0] - h * d2[0] + 2.0 * h * h * d3;
    let c0 = y[0];
    let p = |t: f64| ((c3 * t + c2) * t + c1) * t + c0;
    // split at the critical points
    let mut knots = vec![0.0, 1.0];
    let (a, b, c) = (3.0 * c3, 2.0 * c2, c1);
    if a.abs() > 1e-300 {
        let disc = b * b - 4.0 * a * c;
        if disc >= 0.0 {
            let s = disc.sqrt();
            for r in [(-b - s) / (2.0 * a), (-b + s) / (2.0 * a)] {
                if r > 0.0 && r < 1.0 {
                    knots.push(r);
                }
            }
        }
    } else if b.abs() > 1e-300 {
        let r = -c / b;
        if r > 0.0 && r < 1.0 {
            knots.push(r);
        }
    }
    knots.sort_by(f64::total_cmp);
    let mut roots = Vec::new();
    for w in knots.windows(2) {
        let (mut lo, mut hi) = (w[0], w[1]);
        let (plo, phi) = (p(lo), p(hi));
        if plo == 0.0 {
            roots.push(lo);
            continue;
        }
        if plo.signum() == phi.signum() {
            // near-tangency: keep the knot where |p| is smallest as a candidate
            roots.push(if plo.abs() < phi.abs() { lo } else { hi });
            continue;
        }
        for _ in 0..80 {
            let m = 0.5 * (lo + hi);
            if p(m).signum() == plo.signum() {
                lo = m;
            } else {
                hi = m;
            }
        }
        roots.push(0.5 * (lo + hi));
    }
    if p(1.0) == 0.0 {
        roots.push(1.0);
    }
    roots
}

fn point_triangle_dist(p: V3, a: V3, b: V3, c: V3) -> f64 {
    // Ericson, closest point on triangle
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return norm(ap);
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return norm(bp);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return norm(sub(p, add(a, scl(ab, v))));
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return norm(cp);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return norm(sub(p, add(a, scl(ac, w))));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return norm(sub(p, add(b, scl(sub(c, b), w))));
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    norm(sub(p, add(a, add(scl(ab, v), scl(ac, w)))))
}

fn seg_seg_dist(p1: V3, q1: V3, p2: V3, q2: V3) -> f64 {
    let d1 = sub(q1, p1);
    let d2 = sub(q2, p2);
    let r = sub(p1, p2);
    let a = dot(d1, d1);
    let e = dot(d2, d2);
    let f = dot(d2, r);
    let (s, t);
    if a <= 1e-300 && e <= 1e-300 {
        return norm(r);
    }
    if a <= 1e-300 {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = dot(d1, r);
        if e <= 1e-300 {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = dot(d1, d2);
            let denom = a * e - b * b;
            let mut ss = if denom > 0.0 { ((b * f - c * e) / denom).clamp(0.0, 1.0) } else { 0.0 };
            let mut tt = (b * ss + f) / e;
            if tt < 0.0 {
                tt = 0.0;
                ss = (-c / a).clamp(0.0, 1.0);
            } else if tt > 1.0 {
                tt = 1.0;
                ss = ((b - c) / a).clamp(0.0, 1.0);
            }
            s = ss;
            t = tt;
        }
    }
    norm(sub(add(p1, scl(d1, s)), add(p2, scl(d2, t))))
}

#[derive(Debug, Clone)]
struct Crossing {
    kind: &'static str,
    a: u32,
    b: u32,
    t: f64,
    dist: f64,
    /// side change of the coplanarity function between t = 0 and 1 (odd number of crossings)
    side_change: bool,
}

struct Grid {
    h: f64,
    cells: HashMap<[i64; 3], Vec<u32>>,
}
impl Grid {
    fn key(&self, p: V3) -> [i64; 3] {
        p.map(|c| (c / self.h).floor() as i64)
    }
    fn new(h: f64, boxes: &[(V3, V3)]) -> Grid {
        let mut g = Grid { h, cells: HashMap::new() };
        for (i, (lo, hi)) in boxes.iter().enumerate() {
            let (a, b) = (g.key(*lo), g.key(*hi));
            for x in a[0]..=b[0] {
                for y in a[1]..=b[1] {
                    for z in a[2]..=b[2] {
                        g.cells.entry([x, y, z]).or_default().push(i as u32);
                    }
                }
            }
        }
        g
    }
    fn query(&self, lo: V3, hi: V3, out: &mut Vec<u32>) {
        out.clear();
        let (a, b) = (self.key(lo), self.key(hi));
        for x in a[0]..=b[0] {
            for y in a[1]..=b[1] {
                for z in a[2]..=b[2] {
                    if let Some(v) = self.cells.get(&[x, y, z]) {
                        out.extend_from_slice(v);
                    }
                }
            }
        }
        out.sort_unstable();
        out.dedup();
    }
}

fn bbox(pts: &[V3]) -> (V3, V3) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in pts {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    (lo, hi)
}
fn overlap(a: &(V3, V3), b: &(V3, V3)) -> bool {
    (0..3).all(|k| a.0[k] <= b.1[k] && b.0[k] <= a.1[k])
}

fn edges_of(faces: &[[u32; 3]]) -> Vec<[u32; 2]> {
    let mut e: Vec<[u32; 2]> = faces
        .iter()
        .flat_map(|f| [[f[0], f[1]], [f[1], f[2]], [f[2], f[0]]])
        .map(|[a, b]| if a < b { [a, b] } else { [b, a] })
        .collect();
    e.sort_unstable();
    e.dedup();
    e
}

/// Every vertex-triangle and edge-edge contact of the linear motion `a -> b` (same mesh).
fn swept(faces: &[[u32; 3]], a: &[V3], b: &[V3], l0: f64) -> Vec<Crossing> {
    let tol = 1e-7 * l0;
    let edges = edges_of(faces);
    let tbox: Vec<(V3, V3)> = faces
        .iter()
        .map(|f| bbox(&[a[f[0] as usize], a[f[1] as usize], a[f[2] as usize], b[f[0] as usize], b[f[1] as usize], b[f[2] as usize]]))
        .collect();
    let ebox: Vec<(V3, V3)> = edges
        .iter()
        .map(|e| bbox(&[a[e[0] as usize], a[e[1] as usize], b[e[0] as usize], b[e[1] as usize]]))
        .collect();
    let h = 2.0 * l0;
    let tg = Grid::new(h, &tbox);
    let eg = Grid::new(h, &ebox);
    let mut out: Vec<Crossing> = (0..a.len())
        .into_par_iter()
        .map_init(Vec::new, |buf, v| {
            let vb = bbox(&[a[v], b[v]]);
            tg.query(vb.0, vb.1, buf);
            let mut res = Vec::new();
            for &f in buf.iter() {
                let t = faces[f as usize];
                if t.contains(&(v as u32)) || !overlap(&vb, &tbox[f as usize]) {
                    continue;
                }
                let pos = |s: f64, i: usize| lerp(a[i], b[i], s);
                let g = |s: f64| {
                    let (p, x, y, z) = (pos(s, v), pos(s, t[0] as usize), pos(s, t[1] as usize), pos(s, t[2] as usize));
                    dot(cross(sub(y, x), sub(z, x)), sub(p, x))
                };
                let side_change = g(0.0).signum() != g(1.0).signum();
                for r in cubic_roots(g) {
                    let d = point_triangle_dist(pos(r, v), pos(r, t[0] as usize), pos(r, t[1] as usize), pos(r, t[2] as usize));
                    if d <= tol {
                        res.push(Crossing { kind: "VT", a: v as u32, b: f, t: r, dist: d, side_change });
                    }
                }
            }
            res
        })
        .flatten()
        .collect();
    let ee: Vec<Crossing> = (0..edges.len())
        .into_par_iter()
        .map_init(Vec::new, |buf, i| {
            let e = edges[i];
            eg.query(ebox[i].0, ebox[i].1, buf);
            let mut res = Vec::new();
            for &j in buf.iter() {
                if j as usize <= i {
                    continue;
                }
                let o = edges[j as usize];
                if o.contains(&e[0]) || o.contains(&e[1]) || !overlap(&ebox[i], &ebox[j as usize]) {
                    continue;
                }
                let pos = |s: f64, k: u32| lerp(a[k as usize], b[k as usize], s);
                let g = |s: f64| {
                    let (p, q, r, w) = (pos(s, e[0]), pos(s, e[1]), pos(s, o[0]), pos(s, o[1]));
                    dot(cross(sub(q, p), sub(w, r)), sub(r, p))
                };
                let side_change = g(0.0).signum() != g(1.0).signum();
                for r in cubic_roots(g) {
                    let d = seg_seg_dist(pos(r, e[0]), pos(r, e[1]), pos(r, o[0]), pos(r, o[1]));
                    if d <= tol {
                        res.push(Crossing { kind: "EE", a: i as u32, b: j, t: r, dist: d, side_change });
                    }
                }
            }
            res
        })
        .flatten()
        .collect();
    out.extend(ee);
    out
}

/// Graph distance (rings) between vertex sets, capped: BFS from `src` up to `k` rings.
fn rings_within(adj: &[Vec<u32>], src: &[u32], k: usize) -> Vec<u32> {
    let mut seen: HashMap<u32, usize> = HashMap::new();
    let mut q = VecDeque::new();
    for &s in src {
        seen.insert(s, 0);
        q.push_back(s);
    }
    while let Some(v) = q.pop_front() {
        let d = seen[&v];
        if d == k {
            continue;
        }
        for &w in &adj[v as usize] {
            if !seen.contains_key(&w) {
                seen.insert(w, d + 1);
                q.push_back(w);
            }
        }
    }
    seen.into_keys().collect()
}

/// The smallest distance between two triangles more than `k` rings apart ("sheet distance").
fn sheet_distance(faces: &[[u32; 3]], p: &[V3], l0: f64, k: usize) -> f64 {
    let nv = p.len();
    let mut adj = vec![Vec::new(); nv];
    for e in edges_of(faces) {
        adj[e[0] as usize].push(e[1]);
        adj[e[1] as usize].push(e[0]);
    }
    let r = 3.0 * l0;
    let boxes: Vec<(V3, V3)> = faces
        .iter()
        .map(|f| {
            let (lo, hi) = bbox(&[p[f[0] as usize], p[f[1] as usize], p[f[2] as usize]]);
            (lo.map(|c| c - r / 2.0), hi.map(|c| c + r / 2.0))
        })
        .collect();
    let g = Grid::new(2.0 * l0, &boxes);
    (0..faces.len())
        .into_par_iter()
        .map_init(Vec::new, |buf, i| {
            let f = faces[i];
            g.query(boxes[i].0, boxes[i].1, buf);
            let cand: Vec<u32> = buf.iter().copied().filter(|&j| (j as usize) > i && overlap(&boxes[i], &boxes[j as usize])).collect();
            if cand.is_empty() {
                return f64::INFINITY;
            }
            let near: std::collections::HashSet<u32> = rings_within(&adj, &f, k).into_iter().collect();
            let ti = f.map(|v| p[v as usize]);
            let mut best = f64::INFINITY;
            for j in cand {
                let o = faces[j as usize];
                if o.iter().any(|v| near.contains(v)) {
                    continue;
                }
                let d = ccd::triangle_distance(ti, o.map(|v| p[v as usize]));
                best = best.min(d);
            }
            best
        })
        .reduce(|| f64::INFINITY, f64::min)
}

fn cmd_check(args: &[String]) {
    let dir = PathBuf::from(&args[0]);
    let threads: usize = opt(args, "--threads").unwrap_or(8);
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    let guard_gap: Option<f64> = opt(args, "--guard-gap");
    let mut files: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "bin"))
        .collect();
    files.sort();
    let from: u64 = opt(args, "--from").unwrap_or(0);
    let to: u64 = opt(args, "--to").unwrap_or(u64::MAX);
    let first = read_dump(&files[0]);
    files.retain(|f| {
        let n: u64 = f.file_stem().unwrap().to_str().unwrap()[5..].parse().unwrap();
        n >= from && n <= to
    });
    let l0: f64 = opt(args, "--l0").unwrap_or_else(|| {
        let e = edges_of(&first.faces0);
        e.iter().map(|e| norm(sub(first.p0[e[0] as usize], first.p0[e[1] as usize]))).sum::<f64>() / e.len() as f64
    });
    println!("l0={l0:.6e} steps={}", files.len());
    println!("step,accepted,tau,lib_ccd_tau_max_gap1e-4,max_move_l0,ls_contacts,ls_side_changes,proj_contacts,proj_max_move_l0,sheet_d_l0_proj,remesh_replay_ok,remesh_max_disp_l0,end_intersections,first_ls_contact");
    for f in &files {
        let d = read_dump(f);
        let max_move = d.p0.iter().zip(&d.moved).map(|(a, b)| norm(sub(*a, *b))).fold(0.0, f64::max);
        let ls = swept(&d.faces0, &d.p0, &d.moved, l0);
        let pj = swept(&d.faces0, &d.moved, &d.proj, l0);
        let pj_move = d.moved.iter().zip(&d.proj).map(|(a, b)| norm(sub(*a, *b))).fold(0.0, f64::max);
        let s0 = Surface::new(d.p0.len(), d.faces0.clone()).unwrap();
        let sheet = sheet_distance(&d.faces0, &d.proj, l0, 4);
        // replay the remesh pass (no pins in this scene)
        let mut guard = guard_gap.map(|g| {
            ddg_repulsive_surfaces::ccd::CcdCollapseGuard::new(&CcdOptions { gap: g, ..Default::default() }, l0).unwrap()
        });
        let (out, _) = remesh(
            &s0,
            &d.proj,
            &[],
            l0,
            &RemeshOptions::default(),
            guard.as_mut().map(|g| g as &mut dyn ddg_repulsive_surfaces::remesh::CollapseGuard),
        )
        .unwrap();
        let replay_ok = out.surface.faces() == d.faces1.as_slice() && out.positions == d.p1;
        let mut rmax = 0.0f64;
        for v in 0..out.positions.len() {
            let (src, w) = out.vertex_map.entry(v as u32);
            let mut q = [0.0; 3];
            for (s, wt) in src.iter().zip(w) {
                q = add(q, scl(d.proj[*s as usize], *wt));
            }
            rmax = rmax.max(norm(sub(q, out.positions[v])));
        }
        let s1 = Surface::new(d.p1.len(), d.faces1.clone()).unwrap();
        let endint = ccd::intersecting_pairs(&s1, &d.p1).map_or(usize::MAX, |v| v.len());
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        let vel: Vec<V3> = d.p0.iter().zip(&d.moved).map(|(a, b)| sub(*b, *a)).collect();
        let lib = ccd::first_contact(&pool, &s0, &d.p0, &vel, 1.0, &[], &CcdOptions { gap: 1e-4, reference_length: Some(l0), close_factor: 0.0, ..Default::default() })
            .map_or(String::from("err"), |r| format!("{:.4}", r.tau_max));
        let firstc = ls
            .iter()
            .min_by(|a, b| a.t.total_cmp(&b.t))
            .map_or(String::new(), |c| format!("{} {} {} t={:.4} d={:.2e} side_change={}", c.kind, c.a, c.b, c.t, c.dist, c.side_change));
        println!(
            "{},{},{:.4e},{},{:.3},{},{},{},{:.3e},{:.3},{},{:.3},{},{}",
            d.iteration,
            d.accepted,
            d.tau,
            lib,
            max_move / l0,
            ls.len(),
            ls.iter().filter(|c| c.side_change).count(),
            pj.len(),
            pj_move / l0,
            sheet / l0,
            replay_ok,
            rmax / l0,
            endint,
            firstc
        );
        if !ls.is_empty() {
            let mut vs: Vec<u32> = ls.iter().filter(|c| c.kind == "VT").map(|c| c.a).collect();
            vs.sort_unstable();
            vs.dedup();
            let mut fs: Vec<u32> = ls.iter().filter(|c| c.kind == "VT").map(|c| c.b).collect();
            fs.sort_unstable();
            fs.dedup();
            let cen = |pts: Vec<V3>| { let n = pts.len().max(1) as f64; pts.iter().fold([0.0; 3], |a, p| add(a, *p)).map(|c| c / n) };
            let ts: Vec<f64> = ls.iter().map(|c| c.t).collect();
            println!(
                "  summary step {}: {} VT, {} EE; {} moving vertices through {} faces; t in [{:.3}, {:.3}] of the trial; vertex centroid {:?}, face centroid {:?}; vertex ids {}..{}, face ids {}..{}",
                d.iteration,
                ls.iter().filter(|c| c.kind == "VT").count(),
                ls.iter().filter(|c| c.kind == "EE").count(),
                vs.len(),
                fs.len(),
                ts.iter().copied().fold(f64::INFINITY, f64::min),
                ts.iter().copied().fold(0.0, f64::max),
                cen(vs.iter().map(|&v| d.p0[v as usize]).collect()),
                cen(fs.iter().map(|&f| d.p0[d.faces0[f as usize][0] as usize]).collect()),
                vs.first().unwrap_or(&0), vs.last().unwrap_or(&0), fs.first().unwrap_or(&0), fs.last().unwrap_or(&0)
            );
        }
        if !ls.is_empty() || !pj.is_empty() {
            for c in ls.iter().chain(&pj).take(6) {
                let desc = match c.kind {
                    "VT" => format!("vertex {} at {:?} vs face {:?}", c.a, d.p0[c.a as usize], d.faces0[c.b as usize]),
                    _ => {
                        let e = edges_of(&d.faces0);
                        format!("edge {:?} vs edge {:?}", e[c.a as usize], e[c.b as usize])
                    }
                };
                println!("  crossing step {} {} t={:.5} d={:.2e} side_change={} {}", d.iteration, c.kind, c.t, c.dist, c.side_change, desc);
            }
        }
        std::io::stdout().flush().unwrap();
    }
}

// ---------------------------------------------------------------------------------------------
// topology
// ---------------------------------------------------------------------------------------------

/// Exact linking number of two closed polygons (Klenin & Langowski 2000).
fn linking(a: &[V3], b: &[V3]) -> f64 {
    let mut s = 0.0;
    for i in 0..a.len() {
        let (p1, p2) = (a[i], a[(i + 1) % a.len()]);
        for j in 0..b.len() {
            let (p3, p4) = (b[j], b[(j + 1) % b.len()]);
            let r13 = sub(p3, p1);
            let r14 = sub(p4, p1);
            let r23 = sub(p3, p2);
            let r24 = sub(p4, p2);
            let nn = |u: V3, v: V3| {
                let c = cross(u, v);
                let l = norm(c);
                if l > 0.0 { scl(c, 1.0 / l) } else { [0.0; 3] }
            };
            let n1 = nn(r13, r14);
            let n2 = nn(r14, r24);
            let n3 = nn(r24, r23);
            let n4 = nn(r23, r13);
            let as_ = |x: f64| x.clamp(-1.0, 1.0).asin();
            let om = as_(dot(n1, n2)) + as_(dot(n2, n3)) + as_(dot(n3, n4)) + as_(dot(n4, n1));
            let sg = dot(cross(sub(p4, p3), sub(p2, p1)), r13);
            s += om * sg.signum();
        }
    }
    s / (4.0 * std::f64::consts::PI)
}

/// Tree-cotree generator loops (vertex cycles).
fn generators(faces: &[[u32; 3]], nv: usize) -> Vec<Vec<u32>> {
    let edges = edges_of(faces);
    let eid: HashMap<[u32; 2], usize> = edges.iter().enumerate().map(|(i, e)| (*e, i)).collect();
    let key = |a: u32, b: u32| if a < b { [a, b] } else { [b, a] };
    let mut adj = vec![Vec::new(); nv];
    for e in &edges {
        adj[e[0] as usize].push(e[1]);
        adj[e[1] as usize].push(e[0]);
    }
    // primal BFS tree
    let mut parent = vec![u32::MAX; nv];
    let mut in_tree = vec![false; edges.len()];
    let root = 0u32;
    parent[0] = 0;
    let mut q = VecDeque::from([root]);
    while let Some(v) = q.pop_front() {
        for &w in &adj[v as usize] {
            if parent[w as usize] == u32::MAX {
                parent[w as usize] = v;
                in_tree[eid[&key(v, w)]] = true;
                q.push_back(w);
            }
        }
    }
    // dual tree over the remaining edges
    let mut ef: Vec<Vec<usize>> = vec![Vec::new(); edges.len()];
    for (fi, f) in faces.iter().enumerate() {
        for k in 0..3 {
            ef[eid[&key(f[k], f[(k + 1) % 3])]].push(fi);
        }
    }
    let mut fseen = vec![false; faces.len()];
    let mut in_cotree = vec![false; edges.len()];
    fseen[0] = true;
    let mut q = VecDeque::from([0usize]);
    while let Some(fi) = q.pop_front() {
        let f = faces[fi];
        for k in 0..3 {
            let e = eid[&key(f[k], f[(k + 1) % 3])];
            if in_tree[e] {
                continue;
            }
            for &g in &ef[e] {
                if !fseen[g] {
                    fseen[g] = true;
                    in_cotree[e] = true;
                    q.push_back(g);
                }
            }
        }
    }
    let path = |mut v: u32| {
        let mut p = vec![v];
        while v != root {
            v = parent[v as usize];
            p.push(v);
        }
        p
    };
    let mut loops = Vec::new();
    for (i, e) in edges.iter().enumerate() {
        if in_tree[i] || in_cotree[i] {
            continue;
        }
        let pa = path(e[0]);
        let pb = path(e[1]);
        // trim the common tail
        let (mut ia, mut ib) = (pa.len(), pb.len());
        while ia > 1 && ib > 1 && pa[ia - 2] == pb[ib - 2] {
            ia -= 1;
            ib -= 1;
        }
        // loop: e0 -> ... -> lca -> ... -> e1 -> e0
        let mut l: Vec<u32> = pa[..ia].to_vec();
        l.extend(pb[..ib - 1].iter().rev());
        loops.push(l);
    }
    loops
}

fn vertex_normals(faces: &[[u32; 3]], p: &[V3]) -> Vec<V3> {
    let mut n = vec![[0.0; 3]; p.len()];
    for f in faces {
        let [a, b, c] = f.map(|v| p[v as usize]);
        let fn_ = cross(sub(b, a), sub(c, a));
        for &v in f {
            n[v as usize] = add(n[v as usize], fn_);
        }
    }
    n.iter().map(|x| scl(*x, 1.0 / norm(*x))).collect()
}

fn det(m: &[Vec<f64>]) -> f64 {
    let n = m.len();
    let mut a: Vec<Vec<f64>> = m.to_vec();
    let mut d = 1.0;
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs())).unwrap();
        if a[p][c] == 0.0 {
            return 0.0;
        }
        if p != c {
            a.swap(p, c);
            d = -d;
        }
        d *= a[c][c];
        for r in c + 1..n {
            let f = a[r][c] / a[c][c];
            for k in c..n {
                a[r][k] -= f * a[c][k];
            }
        }
    }
    d
}

fn topo_of(faces: &[[u32; 3]], p: &[V3], label: &str) {
    let nv = p.len();
    let ne = edges_of(faces).len();
    let chi = nv as i64 - ne as i64 + faces.len() as i64;
    let (lo, hi) = bbox(p);
    let mut area = 0.0;
    let mut vol = 0.0;
    for f in faces {
        let [a, b, c] = f.map(|v| p[v as usize]);
        area += 0.5 * norm(cross(sub(b, a), sub(c, a)));
        vol += dot(a, cross(b, c)) / 6.0;
    }
    let e = edges_of(faces);
    let l0 = e.iter().map(|e| norm(sub(p[e[0] as usize], p[e[1] as usize]))).sum::<f64>() / e.len() as f64;
    let s = Surface::new(nv, faces.to_vec()).unwrap();
    let mind = ccd::min_distance(&s, p).unwrap_or(f64::NAN);
    let nint = ccd::intersecting_pairs(&s, p).map_or(usize::MAX, |v| v.len());
    println!(
        "{label}: V={nv} F={} chi={chi} genus={} extents=[{:.3},{:.3},{:.3}] area={area:.4} volume={vol:.4} mean_edge={l0:.4e} min_dist_nonadjacent={mind:.3e} intersecting_pairs={nint}",
        faces.len(),
        (2 - chi) / 2,
        hi[0] - lo[0],
        hi[1] - lo[1],
        hi[2] - lo[2]
    );
    let loops = generators(faces, nv);
    let n = vertex_normals(faces, p);
    let eps = 0.05 * l0.min(if mind.is_finite() { mind } else { l0 });
    let curve = |l: &Vec<u32>, s: f64| -> Vec<V3> { l.iter().map(|&v| add(p[v as usize], scl(n[v as usize], s * eps))).collect() };
    let g = loops.len();
    let mut v = vec![vec![0.0; g]; g];
    let mut maxdev: f64 = 0.0;
    for i in 0..g {
        for j in 0..g {
            let x = linking(&curve(&loops[i], 0.0), &curve(&loops[j], 1.0));
            maxdev = maxdev.max((x - x.round()).abs());
            v[i][j] = x.round();
        }
    }
    println!("  generators={g} loop_lengths={:?} pushoff_eps={eps:.2e} max_non_integer={maxdev:.2e}", loops.iter().map(Vec::len).collect::<Vec<_>>());
    for row in &v {
        println!("  V {:?}", row.iter().map(|x| *x as i64).collect::<Vec<_>>());
    }
    let inter: Vec<Vec<f64>> = (0..g).map(|i| (0..g).map(|j| v[i][j] - v[j][i]).collect()).collect();
    // det(V - t V^T) at t = 0..=g, Newton interpolation -> integer coefficients
    let vals: Vec<f64> = (0..=g)
        .map(|t| {
            let t = t as f64;
            let m: Vec<Vec<f64>> = (0..g).map(|i| (0..g).map(|j| v[i][j] - t * v[j][i]).collect()).collect();
            det(&m)
        })
        .collect();
    // solve Vandermonde
    let k = g + 1;
    let mut a: Vec<Vec<f64>> = (0..k).map(|i| (0..k).map(|j| (i as f64).powi(j as i32)).collect()).collect();
    let mut y = vals.clone();
    for c in 0..k {
        let p = (c..k).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs())).unwrap();
        a.swap(p, c);
        y.swap(p, c);
        for r in 0..k {
            if r != c {
                let f = a[r][c] / a[c][c];
                for q in c..k {
                    a[r][q] -= f * a[c][q];
                }
                y[r] -= f * y[c];
            }
        }
    }
    let coef: Vec<i64> = (0..k).map(|i| (y[i] / a[i][i]).round() as i64).collect();
    println!("  det(V-V^T)={} det(V - t V^T) coefficients (t^0..t^{g}) = {coef:?}", det(&inter).round());
}

fn cmd_topo(args: &[String]) {
    for path in args {
        if path.ends_with(".bin") {
            let d = read_dump(std::path::Path::new(path));
            topo_of(&d.faces0, &d.p0, &format!("{path} (step start)"));
            topo_of(&d.faces1, &d.p1, &format!("{path} (step end)"));
            continue;
        }
        let m = ddg_repulsive_surfaces::read_obj(path).expect("obj");
        let (s, p) = m;
        topo_of(s.faces(), &p, path);
    }
}

// ---------------------------------------------------------------------------------------------
// immersed input (Klein bottle)
// ---------------------------------------------------------------------------------------------

fn ring_dist(adj: &[Vec<u32>], src: &[u32]) -> Vec<u32> {
    let mut d = vec![u32::MAX; adj.len()];
    let mut q = VecDeque::new();
    for &s in src {
        d[s as usize] = 0;
        q.push_back(s);
    }
    while let Some(v) = q.pop_front() {
        for &w in &adj[v as usize] {
            if d[w as usize] == u32::MAX {
                d[w as usize] = d[v as usize] + 1;
                q.push_back(w);
            }
        }
    }
    d
}

/// Immersed input: for every dumped step of a run WITHOUT CCD, the intersection set at t = 0,
/// the library CCD (default gap) on the step's motion, and the exact crossings of the motion
/// classified by their ring distance to the intersection set (k-ring exclusion rule).
fn cmd_klein(args: &[String]) {
    let dir = PathBuf::from(&args[0]);
    let threads: usize = opt(args, "--threads").unwrap_or(8);
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    let mut files: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "bin"))
        .collect();
    files.sort();
    let every: usize = opt(args, "--every").unwrap_or(1);
    let first = read_dump(&files[0]);
    let e = edges_of(&first.faces0);
    let l0 = e.iter().map(|e| norm(sub(first.p0[e[0] as usize], first.p0[e[1] as usize]))).sum::<f64>() / e.len() as f64;
    println!("l0={l0:.5e}");
    println!(
        "step,tau,max_move_l0,int_pairs_t0,int_vertices_t0,lib_tau_max(gap0.05),lib_close,lib_touching,lib_contact,contact_ring,\
         crossings,crossings_ring>0,>1,>2,>3,max_ring_of_crossing"
    );
    for f in files.iter().step_by(every) {
        let d = read_dump(f);
        let s0 = Surface::new(d.p0.len(), d.faces0.clone()).unwrap();
        let ip = ccd::intersecting_pairs(&s0, &d.p0).unwrap();
        let mut iv: Vec<u32> = ip
            .iter()
            .flat_map(|[a, b]| d.faces0[*a as usize].into_iter().chain(d.faces0[*b as usize]))
            .collect();
        iv.sort_unstable();
        iv.dedup();
        let mut adj = vec![Vec::new(); d.p0.len()];
        for e in s0.edges() {
            adj[e[0] as usize].push(e[1]);
            adj[e[1] as usize].push(e[0]);
        }
        let rd = if iv.is_empty() { vec![u32::MAX; d.p0.len()] } else { ring_dist(&adj, &iv) };
        let ring_of = |vs: &[u32]| vs.iter().map(|&v| rd[v as usize]).min().unwrap_or(u32::MAX);
        let vel: Vec<V3> = d.p0.iter().zip(&d.moved).map(|(a, b)| sub(*b, *a)).collect();
        let max_move = vel.iter().map(|v| norm(*v)).fold(0.0, f64::max);
        let rep = ccd::first_contact(
            &pool,
            &s0,
            &d.p0,
            &vel,
            1.0,
            &[],
            &CcdOptions {
                reference_length: Some(l0),
                ..Default::default()
            },
        )
        .unwrap();
        let (cdesc, cring) = match rep.contact.map(|c| c.primitive) {
            Some(ccd::Primitive::VertexTriangle { vertex, triangle }) => {
                let r = ring_of(&[vertex]).max(ring_of(&d.faces0[triangle as usize]));
                (format!("VT v{vertex} f{triangle}"), i64::from(r))
            }
            Some(ccd::Primitive::EdgeEdge { a, b }) => {
                let (ea, eb) = (s0.edges()[a as usize], s0.edges()[b as usize]);
                (format!("EE {ea:?} {eb:?}"), i64::from(ring_of(&ea).max(ring_of(&eb))))
            }
            Some(p) => (format!("{p:?}"), -1),
            None => (String::new(), -1),
        };
        let cr = swept(&d.faces0, &d.p0, &d.moved, l0);
        let edges = edges_of(&d.faces0);
        let rings: Vec<u32> = cr
            .iter()
            .map(|c| match c.kind {
                "VT" => ring_of(&[c.a]).max(ring_of(&d.faces0[c.b as usize])),
                _ => ring_of(&edges[c.a as usize]).max(ring_of(&edges[c.b as usize])),
            })
            .collect();
        let cnt = |k: u32| rings.iter().filter(|&&r| r > k).count();
        println!(
            "{},{:.4e},{:.3},{},{},{:.4},{},{},{},{},{},{},{},{},{},{}",
            d.iteration,
            d.tau,
            max_move / l0,
            ip.len(),
            iv.len(),
            rep.tau_max,
            rep.close_pairs,
            rep.touching_pairs,
            cdesc,
            cring,
            cr.len(),
            cnt(0),
            cnt(1),
            cnt(2),
            cnt(3),
            rings.iter().copied().filter(|&r| r != u32::MAX).max().map_or(-1, i64::from)
        );
        std::io::stdout().flush().unwrap();
    }
}
