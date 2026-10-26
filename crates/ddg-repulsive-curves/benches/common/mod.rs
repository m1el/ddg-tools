//! Shared helpers of the benchmarks (`SPEC.md` §17k): curve loading and refinement, the load
//! average, peak memory, a tiny argument parser. Each bench is a plain `harness = false` main
//! that prints CSV rows.
#![allow(dead_code)]

use std::collections::HashMap;
use std::path::PathBuf;

pub mod components;

use ddg_repulsive_curves::{CurveNetwork, read_curve_obj};

pub type P = [f64; 3];

pub fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/repulsive-curves")
}

/// Freedman–He–Wang with a nominal 600, 1200, 2400 edges (the files have 592, 1184, 2368
/// vertices) as stored, or a uniform refinement (every edge split at its midpoint, repeatedly)
/// of the 2400 curve to the nominal 4800, 9600, 19200 (4736, 9472, 18944 vertices).
pub fn freedman(nominal: usize) -> Result<(CurveNetwork, Vec<P>), String> {
    let (file, doublings) = match nominal {
        600 | 1200 | 2400 => (nominal, 0),
        4800 => (2400, 1),
        9600 => (2400, 2),
        19200 => (2400, 3),
        _ => {
            return Err(format!(
                "no Freedman-He-Wang curve of nominal size {nominal}"
            ));
        }
    };
    let path = data_dir().join(format!("FreedmanHeWang/freedman-he-wang-{file}.obj"));
    let (mut net, mut pos) =
        read_curve_obj(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    for _ in 0..doublings {
        (net, pos) = split_edges(&net, &pos);
    }
    Ok((net, pos))
}

/// Every edge split at its midpoint (new vertices appended).
pub fn split_edges(net: &CurveNetwork, pos: &[P]) -> (CurveNetwork, Vec<P>) {
    let n = net.n_vertices();
    let mut p = pos.to_vec();
    let mut e = Vec::with_capacity(2 * net.edges().len());
    for &[a, b] in net.edges() {
        let m = p.len() as u32;
        let (pa, pb) = (pos[a as usize], pos[b as usize]);
        p.push([0, 1, 2].map(|c| 0.5 * (pa[c] + pb[c])));
        e.push([a, m]);
        e.push([m, b]);
    }
    let _ = n;
    (CurveNetwork::new(p.len(), e).unwrap(), p)
}

/// The closed polyline `pos` (vertices in curve order) resampled to `n` vertices equally spaced
/// in arclength.
pub fn resample_closed(pos: &[P], n: usize) -> Vec<P> {
    let m = pos.len();
    let d =
        |a: P, b: P| ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
    let mut cum = vec![0.0];
    for i in 0..m {
        cum.push(cum[i] + d(pos[i], pos[(i + 1) % m]));
    }
    let total = cum[m];
    let mut out = Vec::with_capacity(n);
    let mut seg = 0;
    for k in 0..n {
        let s = total * k as f64 / n as f64;
        while cum[seg + 1] < s {
            seg += 1;
        }
        let t = (s - cum[seg]) / (cum[seg + 1] - cum[seg]);
        let (a, b) = (pos[seg], pos[(seg + 1) % m]);
        out.push([0, 1, 2].map(|c| a[c] + t * (b[c] - a[c])));
    }
    out
}

pub fn closed_network(n: usize) -> CurveNetwork {
    CurveNetwork::new(n, (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect()).unwrap()
}

/// Knot 96 (Knot128 perturbed, 2106 vertices) resampled by arclength to `n` vertices.
pub fn knot96(n: usize) -> Result<(CurveNetwork, Vec<P>), String> {
    let path = data_dir().join("KnotDatasets/Knot128/perturbed/knot0096.obj");
    let (net, pos) = read_curve_obj(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    // vertices of a closed single curve may not be in curve order: walk the edges
    let order = walk_closed(&net)?;
    let ordered: Vec<P> = order.iter().map(|&v| pos[v as usize]).collect();
    let p = resample_closed(&ordered, n);
    Ok((closed_network(n), p))
}

/// The Freedman–He–Wang 2400 curve resampled to `n` vertices (the paper's own sizes, e.g. 1000).
pub fn freedman_resampled(n: usize) -> Result<(CurveNetwork, Vec<P>), String> {
    let (net, pos) = freedman(2400)?;
    let order = walk_closed(&net)?;
    let ordered: Vec<P> = order.iter().map(|&v| pos[v as usize]).collect();
    Ok((closed_network(n), resample_closed(&ordered, n)))
}

pub fn walk_closed(net: &CurveNetwork) -> Result<Vec<u32>, String> {
    let n = net.n_vertices();
    let mut adj = vec![[u32::MAX; 2]; n];
    for &[a, b] in net.edges() {
        for (x, y) in [(a, b), (b, a)] {
            let s = &mut adj[x as usize];
            if s[0] == u32::MAX {
                s[0] = y
            } else if s[1] == u32::MAX {
                s[1] = y
            } else {
                return Err("vertex of degree > 2".into());
            }
        }
    }
    let mut order = vec![0u32];
    let (mut prev, mut cur) = (u32::MAX, 0u32);
    loop {
        let [a, b] = adj[cur as usize];
        let next = if a != prev { a } else { b };
        if next == 0 || order.len() == n {
            break;
        }
        order.push(next);
        prev = cur;
        cur = next;
    }
    if order.len() != n {
        return Err(format!("not one closed curve ({} of {n})", order.len()));
    }
    Ok(order)
}

/// A curve by name: `fhw` (Freedman–He–Wang, sizes 600..19200), `fhw-r` (resampled), `knot`.
pub fn curve(name: &str, n: usize) -> Result<(CurveNetwork, Vec<P>), String> {
    match name {
        "fhw" => freedman(n),
        "fhw-r" => freedman_resampled(n),
        "knot" => knot96(n),
        _ => Err(format!("unknown curve {name}")),
    }
}

pub fn load_avg() -> f64 {
    std::fs::read_to_string("/proc/loadavg")
        .ok()
        .and_then(|s| s.split_whitespace().next().and_then(|x| x.parse().ok()))
        .unwrap_or(f64::NAN)
}

/// Peak resident set size in GB (`VmHWM`).
pub fn peak_rss_gb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM"))
                .and_then(|l| l.split_whitespace().nth(1)?.parse::<f64>().ok())
        })
        .map(|kb| kb / 1.0e6)
        .unwrap_or(f64::NAN)
}

/// `--key value` pairs (cargo's own `--bench` flag is ignored).
pub fn args() -> HashMap<String, String> {
    let v: Vec<String> = std::env::args().skip(1).collect();
    let mut m = HashMap::new();
    let mut i = 0;
    while i < v.len() {
        if let Some(k) = v[i].strip_prefix("--") {
            if k == "bench" {
                i += 1;
                continue;
            }
            let val = v.get(i + 1).cloned().unwrap_or_default();
            m.insert(k.to_string(), val);
            i += 2;
        } else {
            i += 1;
        }
    }
    m
}

pub fn get<T: std::str::FromStr>(a: &HashMap<String, String>, k: &str, d: T) -> T {
    a.get(k).and_then(|s| s.parse().ok()).unwrap_or(d)
}

pub fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if s.is_empty() {
        f64::NAN
    } else {
        s[s.len() / 2]
    }
}
