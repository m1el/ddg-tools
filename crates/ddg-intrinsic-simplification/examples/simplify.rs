//! Coarsens a mesh with the ICE metric (paper Alg. 1) and prints the report; optionally
//! compares the removal order, the costs and the coarse lengths with a run of the reference
//! driver (`runs_tmp/iem/reference-driver.md`). Single thread.
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-simplification --example simplify -- MESH.obj \
//!     [--target N | --ratio r | --min] [--compat] [--w W] [--memoryless] [--one-ring] \
//!     [--no-delaunay] [--out DIR] [--ref DIR] [--repeat K] [--certify]
//! ```
//!
//! * `--compat`: `SimplifyOptions::reference_compat(W, …)` (SPEC §12.2); otherwise the
//!   defaults with `mixture(W)` if `--w` is given.
//! * `--out DIR`: writes `removal.tsv` (`step vertex type key cost`, types as the reference:
//!   0 interior, 1 ear, 2 boundary), `coarse.txt` (per face: three fine ids and three lengths)
//!   and `summary.txt`.
//! * `--ref DIR`: a reference driver output directory: compares the removal sequence step by
//!   step (vertex and cost) and, if the sequences agree, the coarse lengths.
//! * `--repeat K`: runs `K` times and reports the fastest (timings).

use std::collections::HashMap;
use std::fmt::Write as _;
use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_intrinsic_simplification::{
    ChannelWeights, CoarseMesh, CoarseMeshOptions, RemovalKind, Requeue, SimplifyOptions,
    SimplifyReport, Target, simplify,
};

fn arg<T: std::str::FromStr>(args: &[String], name: &str) -> Option<T> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|x| x.parse().ok())
}

fn die(msg: String) -> ! {
    eprintln!("{msg}");
    std::process::exit(1);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(path) = args.get(1) else {
        die("usage: simplify MESH.obj [--target N | --ratio r | --min] [--compat] ...".into());
    };
    let flag = |name: &str| args.iter().any(|a| a == name);
    let target = if let Some(n) = arg::<usize>(&args, "--target") {
        Target::Vertices(n)
    } else if let Some(r) = arg::<f64>(&args, "--ratio") {
        Target::Ratio(r)
    } else if flag("--min") {
        Target::AllRemovable
    } else {
        Target::Vertices(500)
    };
    let w: Option<f64> = arg(&args, "--w");
    let mut options = if flag("--compat") {
        SimplifyOptions::reference_compat(w.unwrap_or(0.0), target)
    } else {
        SimplifyOptions {
            target,
            channels: w.map_or_else(ChannelWeights::curvature, ChannelWeights::mixture),
            ..SimplifyOptions::default()
        }
    };
    if flag("--memoryless") {
        options.metric.memory = false;
    }
    if flag("--one-ring") {
        options.requeue = Requeue::OneRing;
    }
    if flag("--no-delaunay") {
        options.initial_delaunay = false;
    }
    let out: Option<String> = arg(&args, "--out");
    let reference: Option<String> = arg(&args, "--ref");
    options.record_sequence = out.is_some() || reference.is_some();
    let repeat: usize = arg(&args, "--repeat").unwrap_or(1);

    let obj = read_obj(path).unwrap_or_else(|e| die(format!("{path}: {e}")));
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let mut best: Option<(f64, CoarseMesh, SimplifyReport)> = None;
    for _ in 0..repeat.max(1) {
        let mut mesh =
            CoarseMesh::from_positions(&obj.positions, &faces, &[], CoarseMeshOptions::default())
                .unwrap_or_else(|e| die(format!("{path}: {e}")));
        let t = Instant::now();
        let report = simplify(&mut mesh, &options).unwrap_or_else(|e| die(format!("{path}: {e}")));
        let secs = t.elapsed().as_secs_f64();
        if best.as_ref().is_none_or(|b| secs < b.0) {
            best = Some((secs, mesh, report));
        }
    }
    let (secs, mesh, report) = best.unwrap_or_else(|| die("no run".into()));
    let n_removed = report.n_interior + report.n_boundary + report.n_ear;
    let mut summary = String::new();
    let _ = writeln!(
        summary,
        "{path}: {} -> {} vertices ({:?}); removed {n_removed} (interior {}, boundary {}, ear {})",
        report.n_vertices_before,
        report.n_vertices_after,
        report.stop,
        report.n_interior,
        report.n_boundary,
        report.n_ear
    );
    let t = &report.timings;
    let _ = writeln!(
        summary,
        "time {secs:.4} s ({:.0} removals/s; best of {repeat}): initial Delaunay {:.4}, \
         channels {:.4}, initial costs {:.4}, removals {:.4}, requeue {:.4}, snapshots {:.4}",
        n_removed as f64 / secs,
        t.initial_delaunay,
        t.channels,
        t.initial_costs,
        t.removals,
        t.requeue,
        t.snapshots
    );
    let _ = writeln!(
        summary,
        "pops {} (stale {}), evaluations {}, commit failures {:?}, infinite evaluations {:?}, \
         skips at the end {:?}",
        report.n_pops,
        report.n_stale_pops,
        report.n_evaluations,
        report.commit_failures,
        report.infinite_evaluations,
        report.skips
    );
    let _ = writeln!(
        summary,
        "flips: initial Delaunay {}, degree {}, fallback {}, Delaunay {} (refused {}); \
         flattening iterations max {} mean {:.2}; clamped {}",
        report.n_initial_delaunay_flips,
        report.n_degree_flips,
        report.n_fallback_flips,
        report.n_delaunay_flips,
        report.n_delaunay_refused,
        report.flatten_iterations_max,
        report.flatten_iterations_mean,
        report.n_clamped
    );
    let _ = writeln!(summary, "keys (S5): {:?}", report.keys);
    let _ = writeln!(summary, "components: {:?}", report.components);
    // Quality: minimum corner angle, curvature distribution.
    let mut min_angle = f64::INFINITY;
    for h in 0..mesh.complex().n_halfedges() {
        min_angle = min_angle.min(mesh.corner_angle(h));
    }
    let mut curv: Vec<f64> = (0..mesh.n_vertices())
        .map(|v| mesh.curvature(v).abs())
        .collect();
    curv.sort_by(f64::total_cmp);
    let q = |x: f64| curv[((curv.len() - 1) as f64 * x).round() as usize];
    let _ = writeln!(
        summary,
        "min corner angle {:.3} deg; |K| quantiles 50% {:.3e} 90% {:.3e} max {:.3e}; sum |K| {:.6}; \
         total area {:.6}",
        min_angle.to_degrees(),
        q(0.5),
        q(0.9),
        q(1.0),
        curv.iter().sum::<f64>(),
        mesh.total_area()
    );
    match mesh.certify(options.initial_delaunay || options.reference.is_some()) {
        Ok(c) => {
            let _ = writeln!(summary, "certify: ok ({} not Delaunay)", c.n_not_delaunay);
        }
        Err(e) => {
            let _ = writeln!(summary, "certify: FAILED: {e}");
        }
    }
    print!("{summary}");
    let kind = |k: RemovalKind| match k {
        RemovalKind::Interior => 0,
        RemovalKind::Ear => 1,
        RemovalKind::Boundary => 2,
    };
    if let Some(dir) = &out {
        std::fs::create_dir_all(dir).unwrap_or_else(|e| die(format!("{dir}: {e}")));
        let mut s = String::from("step\tvertex\ttype\tkey\tcost\n");
        for (i, r) in report.sequence.iter().enumerate() {
            let _ = writeln!(
                s,
                "{}\t{}\t{}\t{:e}\t{:e}",
                i + 1,
                r.vertex,
                kind(r.kind),
                r.key,
                r.cost
            );
        }
        write(&format!("{dir}/removal.tsv"), &s);
        let mut s = String::new();
        let c = mesh.complex();
        for (f, face) in mesh.faces().iter().enumerate() {
            let l = [0, 1, 2].map(|k| mesh.edge_lengths()[c.edge(3 * f + k)]);
            let _ = writeln!(
                s,
                "{} {} {} {:e} {:e} {:e}",
                face[0], face[1], face[2], l[0], l[1], l[2]
            );
        }
        write(&format!("{dir}/coarse.txt"), &s);
        write(&format!("{dir}/summary.txt"), &summary);
    }
    if let Some(dir) = &reference {
        compare(dir, &report, &mesh, kind);
        if flag("--replay") {
            replay(dir, &obj.positions, &faces, w.unwrap_or(0.0), kind);
        }
    }
}

/// Removes the vertices in the reference's order with the compat metric (public API only)
/// and compares every commit cost and the final lengths: a check of the metric and the
/// transport that does not depend on how near-ties are broken.
fn replay(
    dir: &str,
    positions: &[[f64; 3]],
    faces: &[[u32; 3]],
    w: f64,
    kind: fn(RemovalKind) -> i32,
) {
    use ddg_intrinsic::DelaunayCriterion;
    use ddg_intrinsic_simplification::ChannelState;
    let mut compat = SimplifyOptions::reference_compat(w, Target::AllRemovable);
    if std::env::args().any(|a| a == "--flat-threshold") {
        // Against a reference built with IEM_FLAT_THRESHOLD (equal shares iff |K_i| <= 1e-12).
        compat.metric.flat = ddg_intrinsic_simplification::FlatRule::Threshold(1e-12);
    }
    let r = compat.reference.unwrap_or_default();
    let mesh_options = CoarseMeshOptions {
        delaunay_after_removal: false,
        delaunay: DelaunayCriterion::AngleSum {
            tolerance: r.delaunay_tolerance.unwrap_or(1e-10),
        },
        flatten_tolerance: r.flatten_tolerance.unwrap_or(1e-12),
        ..CoarseMeshOptions::default()
    };
    let mut mesh = CoarseMesh::from_positions(positions, faces, &[], mesh_options)
        .unwrap_or_else(|e| die(format!("{e}")));
    let mut state = ChannelState::new(&mesh, &compat.channels, r.fold_weights)
        .unwrap_or_else(|e| die(format!("{e}")));
    let text = read(&format!("{dir}/removal.tsv"));
    let rows: Vec<(u32, i32, f64)> = text
        .lines()
        .skip(1)
        .filter_map(|l| {
            let x: Vec<&str> = l.split('\t').collect();
            Some((
                x.get(1)?.parse().ok()?,
                x.get(2)?.parse().ok()?,
                x.get(4)?.parse().ok()?,
            ))
        })
        .collect();
    let mut max_rel: f64 = 0.0;
    let mut worst = 0;
    let mut n_over = [0usize; 3];
    let mut rels: Vec<f64> = Vec::with_capacity(rows.len());
    for (step, &(v, t, ref_cost)) in rows.iter().enumerate() {
        let prepared = match mesh.begin_removal(v) {
            Ok(p) => p,
            Err(e) => {
                println!("replay: step {}: {e}; stopping", step + 1);
                return;
            }
        };
        if kind(prepared.kind) != t {
            println!(
                "replay: step {}: vertex {v} is of type {} here, {t} in the reference",
                step + 1,
                kind(prepared.kind)
            );
        }
        let cost = state.apply_removal(&prepared, &compat.metric, r.commit_frames);
        let removal = mesh
            .finish_removal()
            .unwrap_or_else(|e| die(format!("replay: step {}: {e}", step + 1)));
        // R21: every edge of the faces within two edge steps of the merged face.
        let c = mesh.complex();
        let mut ring = vec![removal.face as usize];
        let mut start = 0;
        for _ in 0..2 {
            let end = ring.len();
            for k in start..end {
                for s in 0..3 {
                    if let Some(tw) = c.twin(3 * ring[k] + s)
                        && !ring.contains(&(tw / 3))
                    {
                        ring.push(tw / 3);
                    }
                }
            }
            start = end;
        }
        let mut seeds: Vec<u32> = ring
            .iter()
            .flat_map(|&g| (0..3).map(move |s| c.edge(3 * g + s) as u32))
            .collect();
        seeds.sort_unstable();
        seeds.dedup();
        mesh.flip_to_delaunay_from(&seeds)
            .unwrap_or_else(|e| die(format!("replay: step {}: {e}", step + 1)));
        let rel = (cost - ref_cost).abs() / ref_cost.abs().max(1e-300);
        rels.push(rel);
        for (k, tol) in [1e-12, 1e-9, 1e-6].iter().enumerate() {
            if rel > *tol {
                n_over[k] += 1;
            }
        }
        if rel > 1e-6 && n_over[2] <= 3 {
            println!(
                "  replay step {}: vertex {v} type {t}: ours {cost:e}, reference {ref_cost:e} \
                 ({} neighbours, curvature {:e})",
                step + 1,
                prepared.neighbours.len(),
                prepared.curvature
            );
        }
        if rel > max_rel {
            max_rel = rel;
            worst = step + 1;
        }
    }
    println!(
        "replay of the reference order ({} removals): commit cost relative difference max \
         {max_rel:.3e} (step {worst}); steps above 1e-12: {}, 1e-9: {}, 1e-6: {}",
        rows.len(),
        n_over[0],
        n_over[1],
        n_over[2]
    );
    rels.sort_by(f64::total_cmp);
    let q = |x: f64| {
        rels.get(((rels.len().max(1) - 1) as f64 * x) as usize)
            .copied()
            .unwrap_or(0.0)
    };
    println!(
        "  relative differences: median {:.2e}, 90% {:.2e}, 99% {:.2e}; above 1e-3: {}, above \
         1e-2: {}",
        q(0.5),
        q(0.9),
        q(0.99),
        rels.iter().filter(|&&r| r > 1e-3).count(),
        rels.iter().filter(|&&r| r > 1e-2).count()
    );
    compare_lengths(dir, &mesh);
}

fn write(path: &str, s: &str) {
    std::fs::write(path, s).unwrap_or_else(|e| die(format!("{path}: {e}")));
}

fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|e| die(format!("{path}: {e}")))
}

/// Compares with a reference driver run.
fn compare(dir: &str, report: &SimplifyReport, mesh: &CoarseMesh, kind: fn(RemovalKind) -> i32) {
    let text = read(&format!("{dir}/removal.tsv"));
    let rows: Vec<(u32, i32, f64, f64)> = text
        .lines()
        .skip(1)
        .filter_map(|l| {
            let x: Vec<&str> = l.split('\t').collect();
            Some((
                x.get(1)?.parse().ok()?,
                x.get(2)?.parse().ok()?,
                x.get(3)?.parse().ok()?,
                x.get(4)?.parse().ok()?,
            ))
        })
        .collect();
    let ours = &report.sequence;
    let n = rows.len().min(ours.len());
    let mut first_divergence = None;
    let mut max_key_rel: f64 = 0.0;
    let mut max_cost_rel: f64 = 0.0;
    for k in 0..n {
        let (v, t, key, cost) = rows[k];
        let r = &ours[k];
        if r.vertex != v || kind(r.kind) != t {
            first_divergence = Some(k);
            break;
        }
        let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
        max_key_rel = max_key_rel.max(rel(r.key, key));
        max_cost_rel = max_cost_rel.max(rel(r.cost, cost));
    }
    let agree = first_divergence.unwrap_or(n);
    println!(
        "reference {dir}: {} removals (ours {}); the first {agree} agree (vertex and type); \
         max relative difference over them: key {max_key_rel:.3e}, commit cost {max_cost_rel:.3e}",
        rows.len(),
        ours.len()
    );
    // Same set of removed vertices overall.
    let theirs: std::collections::HashSet<u32> = rows.iter().map(|r| r.0).collect();
    let mine: std::collections::HashSet<u32> = ours.iter().map(|r| r.vertex).collect();
    println!(
        "removed sets: common {}, only reference {}, only ours {}",
        theirs.intersection(&mine).count(),
        theirs.difference(&mine).count(),
        mine.difference(&theirs).count()
    );
    // Agreement up to the order of near-ties: D(k) = vertices among the first k removals of
    // one side and not of the other.
    {
        let mut in_ref = HashMap::new();
        let mut in_ours = HashMap::new();
        let mut d = 0i64;
        let (mut equal_steps, mut last_equal, mut max_d) = (0usize, 0usize, 0i64);
        let mut first_large = None;
        for k in 0..n {
            let a = rows[k].0;
            let b = ours[k].vertex;
            if a != b {
                // a enters the reference prefix; b ours.
                if in_ours.remove(&a).is_some() {
                    d -= 1;
                } else {
                    in_ref.insert(a, ());
                    d += 1;
                }
                if in_ref.remove(&b).is_some() {
                    d -= 1;
                } else {
                    in_ours.insert(b, ());
                    d += 1;
                }
            }
            let dk = d / 2;
            max_d = max_d.max(dk);
            if dk == 0 {
                equal_steps += 1;
                last_equal = k + 1;
            }
            if dk > 10 && first_large.is_none() {
                first_large = Some(k + 1);
            }
        }
        println!(
            "prefix sets: equal after {equal_steps} of {n} steps, last equal at step {last_equal}; \
             largest difference {max_d} vertices; first above 10 at step {first_large:?}; \
             difference at the end {}",
            d / 2
        );
    }
    if let Some(k) = first_divergence {
        let (v, t, key, cost) = rows[k];
        let r = &ours[k];
        println!(
            "first divergence at step {}: reference removes {v} (type {t}, key {key:e}, cost \
             {cost:e}); ours {} (type {}, key {:e}, cost {:e})",
            k + 1,
            r.vertex,
            kind(r.kind),
            r.key,
            r.cost
        );
        // Where does each side remove the other's vertex?
        let pos_ours = ours.iter().position(|x| x.vertex == v);
        let pos_ref = rows.iter().position(|x| x.0 == r.vertex);
        println!(
            "  ours removes {v} at step {:?}; the reference removes {} at step {:?}",
            pos_ours.map(|p| p + 1),
            r.vertex,
            pos_ref.map(|p| p + 1)
        );
        for j in k.saturating_sub(2)..(k + 3).min(n) {
            println!(
                "  step {}: ref {} {:e} | ours {} {:e}",
                j + 1,
                rows[j].0,
                rows[j].2,
                ours[j].vertex,
                ours[j].key
            );
        }
        return;
    }
    compare_lengths(dir, mesh);
}

/// Compares the coarse lengths with the reference's: edges keyed by their fine end points.
fn compare_lengths(dir: &str, mesh: &CoarseMesh) {
    let vidx: Vec<u32> = read(&format!("{dir}/coarse_vIdx.txt"))
        .split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect();
    let f: Vec<u32> = read(&format!("{dir}/coarse_F.txt"))
        .split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect();
    let l: Vec<f64> = read(&format!("{dir}/coarse_l.txt"))
        .split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect();
    let mut theirs: HashMap<(u32, u32), Vec<f64>> = HashMap::new();
    for (face, lens) in f.chunks(3).zip(l.chunks(3)) {
        for s in 0..3 {
            let (a, b) = (vidx[face[s] as usize], vidx[face[(s + 1) % 3] as usize]);
            theirs
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push(lens[s]);
        }
    }
    let mut mine: HashMap<(u32, u32), Vec<f64>> = HashMap::new();
    let c = mesh.complex();
    for (fi, face) in mesh.faces().iter().enumerate() {
        for s in 0..3 {
            let (a, b) = (face[s], face[(s + 1) % 3]);
            mine.entry((a.min(b), a.max(b)))
                .or_default()
                .push(mesh.edge_lengths()[c.edge(3 * fi + s)]);
        }
    }
    let mut max_rel: f64 = 0.0;
    let mut missing = 0;
    for (key, ls) in &theirs {
        match mine.get(key) {
            Some(ms) if ms.len() == ls.len() => {
                let mut a = ls.clone();
                let mut b = ms.clone();
                a.sort_by(f64::total_cmp);
                b.sort_by(f64::total_cmp);
                for (x, y) in a.iter().zip(&b) {
                    max_rel = max_rel.max((x - y).abs() / x);
                }
            }
            _ => missing += 1,
        }
    }
    println!(
        "coarse edges: reference {}, ours {}, edges not matched {missing}; max relative length \
         difference over matched {max_rel:.3e}",
        theirs.len(),
        mine.len()
    );
}
