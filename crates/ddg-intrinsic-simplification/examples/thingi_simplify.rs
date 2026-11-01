//! WP9b, SPEC §10.2: the Thingi10k benchmark of intrinsic simplification (paper Sec. 8.2, Fig. 16-17)
//! for one mesh.
//!
//! ```text
//! thingi_simplify MESH.obj [--id ID] [--no-refine] [--signposts] [--target-ratio R] [--no-mollify]
//! ```
//!
//! The edge lengths are mollified by default (`eps = 1e-5` of the mean edge length, as in the
//! integer-coordinates evaluation): Thingi10k has faces of zero area, which `Geometry` refuses.
//!
//! One process, one mesh; prints one JSON line on stdout (also on failure). With pre-refinement
//! (default: integer coordinates, Chew's rules, 25 degrees; `--signposts`: the signpost
//! representation) or, with `--no-refine`, on the input triangulation. Coarsens as far as
//! possible (`Target::AllRemovable`), default ICE metric; `--target-ratio R` uses
//! `Target::PerComponent(R)` instead.
//!
//! The ratio: with `c` components, `|V|` input vertices and `|V̄|` vertices at the end,
//! `ratio = |V̄| / |V|` (the paper normalizes both counts by `c`, `(|V̄|/c)/(|V|/c)`, which is
//! this). `ratio_mean_component` is the mean over components of `|V̄_k| / |V_k|` (each component
//! against its own input vertices, refinement not counted). The driver `scripts/iem-thingi.sh` runs it in parallel.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::time::Instant;

use ddg_core::io::parse_obj;
use ddg_core::mollify::mollify_edge_lengths;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::PAPER_MOLLIFY_EPSILON;
use ddg_intrinsic_simplification::refine::{
    RefinementOptions, Representation, refine_then_simplify,
};
use ddg_intrinsic_simplification::{
    CoarseMesh, CoarseMeshOptions, SimplifyOptions, SimplifyReport, Target, simplify,
};

#[derive(Default)]
struct Json(String);

impl Json {
    fn sep(&mut self) {
        if !self.0.is_empty() {
            self.0.push(',');
        }
    }
    fn raw(&mut self, k: &str, v: &str) -> &mut Self {
        self.sep();
        let _ = write!(self.0, "\"{k}\":{v}");
        self
    }
    fn num(&mut self, k: &str, v: impl std::fmt::Display) -> &mut Self {
        let s = v.to_string();
        let s = if s == "inf" || s == "-inf" || s == "NaN" {
            "null".to_string()
        } else {
            s
        };
        self.raw(k, &s)
    }
    fn str(&mut self, k: &str, v: &str) -> &mut Self {
        let esc: String = v
            .chars()
            .flat_map(|c| match c {
                '"' => vec!['\\', '"'],
                '\\' => vec!['\\', '\\'],
                '\n' | '\t' | '\r' => vec![' '],
                c => vec![c],
            })
            .take(600)
            .collect();
        self.raw(k, &format!("\"{esc}\""))
    }
    fn map(&mut self, k: &str, m: &BTreeMap<&'static str, usize>) -> &mut Self {
        let body: Vec<String> = m.iter().map(|(a, b)| format!("\"{a}\":{b}")).collect();
        self.raw(k, &format!("{{{}}}", body.join(",")))
    }
}

fn peak_rss_kb() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1).and_then(|x| x.parse().ok()))
        })
        .unwrap_or(0)
}

/// Per component (ordered by smallest vertex id) the number of input vertices.
fn component_sizes(n: usize, faces: &[[u32; 3]]) -> Vec<usize> {
    fn find(p: &mut [u32], mut x: u32) -> u32 {
        while p[x as usize] != x {
            p[x as usize] = p[p[x as usize] as usize];
            x = p[x as usize];
        }
        x
    }
    let mut p: Vec<u32> = (0..n as u32).collect();
    for f in faces {
        let a = find(&mut p, f[0]);
        for &v in &f[1..] {
            let b = find(&mut p, v);
            if a != b {
                let a = find(&mut p, f[0]);
                p[b as usize] = a;
            }
        }
    }
    let mut index = vec![u32::MAX; n];
    let mut sizes = Vec::new();
    for v in 0..n as u32 {
        let r = find(&mut p, v) as usize;
        if index[r] == u32::MAX {
            index[r] = sizes.len() as u32;
            sizes.push(0);
        }
        sizes[index[r] as usize] += 1;
    }
    sizes
}

fn report_json(j: &mut Json, r: &SimplifyReport, comp_in: &[usize], v0: usize) {
    let after_total = r.n_vertices_after;
    j.str("stop", &format!("{:?}", r.stop))
        .num("v_before_simplify", r.n_vertices_before)
        .num("v_final", after_total)
        .num("ratio", after_total as f64 / v0.max(1) as f64);
    if r.components.len() == comp_in.len() && !comp_in.is_empty() {
        let mean = r
            .components
            .iter()
            .zip(comp_in)
            .map(|(c, &n)| c.after as f64 / n.max(1) as f64)
            .sum::<f64>()
            / comp_in.len() as f64;
        j.num("ratio_mean_component", mean);
        let after: Vec<String> = r.components.iter().map(|c| c.after.to_string()).collect();
        // Per component final counts, only for a few components (the long tail is in the sum).
        if after.len() <= 64 {
            j.raw("v_final_per_component", &format!("[{}]", after.join(",")));
        }
    }
    j.num("removed_interior", r.n_interior)
        .num("removed_boundary", r.n_boundary)
        .num("removed_ear", r.n_ear)
        .map("skips", &r.skips)
        .num("skips_total", r.skips.values().sum::<usize>())
        .map("infinite_evaluations", &r.infinite_evaluations)
        .map("commit_failures", &r.commit_failures)
        .num("commit_failures_total", r.commit_failures.values().sum::<usize>())
        .num(
            "flatten_failures",
            r.commit_failures
                .iter()
                .filter(|(k, _)| k.contains("flatten"))
                .map(|(_, v)| *v)
                .sum::<usize>(),
        )
        .num("n_clamped", r.n_clamped)
        .num("n_pops", r.n_pops)
        .num("n_evaluations", r.n_evaluations)
        .num("n_fallback_flips", r.n_fallback_flips)
        .num("n_degree_flips", r.n_degree_flips)
        .num("n_delaunay_flips", r.n_delaunay_flips)
        .num("n_delaunay_refused", r.n_delaunay_refused)
        .num("flatten_iterations_max", r.flatten_iterations_max)
        .num("key_mismatch_clean", r.keys.n_clean_mismatch)
        .num("t_initial_delaunay", r.timings.initial_delaunay)
        .num("t_channels", r.timings.channels)
        .num("t_initial_costs", r.timings.initial_costs)
        .num("t_removals", r.timings.removals)
        .num("t_requeue", r.timings.requeue)
        .num("t_simplify", r.timings.total);
    let removed = r.n_interior + r.n_boundary + r.n_ear;
    j.num("removals", removed);
    if removed > 0 {
        j.num("sec_per_removal", r.timings.total / removed as f64);
    }
}

fn run(path: &str, args: &[String], j: &mut Json) -> Result<(), String> {
    let refine = !args.iter().any(|a| a == "--no-refine");
    let signposts = args.iter().any(|a| a == "--signposts");
    let target = args
        .iter()
        .position(|a| a == "--target-ratio")
        .and_then(|i| args.get(i + 1))
        .and_then(|x| x.parse::<f64>().ok())
        .map(Target::PerComponent)
        .unwrap_or(Target::AllRemovable);
    let t = Instant::now();
    let text = std::fs::read(path).map_err(|e| format!("read: {e}"))?;
    let obj = parse_obj(&String::from_utf8_lossy(&text)).map_err(|e| format!("parse: {e}"))?;
    let faces = obj.triangles().ok_or("not all triangles")?;
    let mesh = TriMesh::new(obj.positions.len(), faces.clone()).map_err(|e| format!("mesh: {e}"))?;
    // Intrinsic mollification (paper Sec. 3.9 of Gillespie et al. 2021, eps = 1e-5 mean edge
    // length), as the integer-coordinates evaluation: Thingi10k has zero-area faces.
    let lengths: Vec<f64> = mesh
        .edges()
        .iter()
        .map(|&[a, b]| {
            let (a, b) = (obj.positions[a as usize], obj.positions[b as usize]);
            ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
        })
        .collect();
    let geometry = if args.iter().any(|a| a == "--no-mollify") {
        Geometry::from_edge_lengths(&mesh, lengths)
    } else {
        let m = mollify_edge_lengths(&mesh, &lengths, PAPER_MOLLIFY_EPSILON)
            .map_err(|e| format!("mollify: {e}"))?;
        j.num("mollify_epsilon", m.epsilon)
            .num("mollify_faces_below", m.n_faces_below_margin);
        Geometry::from_edge_lengths(&mesh, m.edge_lengths)
    }
    .map_err(|e| format!("geometry: {e}"))?;
    drop(text);
    let v0 = mesh.n_vertices();
    let comp_in = component_sizes(v0, &faces);
    j.num("v0", v0)
        .num("f0", mesh.n_faces())
        .num("boundary_edges", mesh.n_boundary_edges())
        .num("components", mesh.n_components())
        .num("chi", mesh.euler_characteristic())
        .str("mode", if !refine { "no-refine" } else if signposts { "refine-signposts" } else { "refine" })
        .num("t_load", t.elapsed().as_secs_f64());
    let so = SimplifyOptions {
        target,
        ..SimplifyOptions::default()
    };
    if refine {
        let ro = RefinementOptions {
            representation: if signposts {
                Representation::Signposts
            } else {
                Representation::DEFAULT
            },
            ..RefinementOptions::default()
        };
        let (coarse, rep) =
            refine_then_simplify(&mesh, &geometry, &[], &ro, &so).map_err(|e| format!("{e}"))?;
        let r = &rep.refinement;
        j.num("refine_added", r.n_added)
            .num("refine_min_angle_deg", r.min_angle.to_degrees())
            .num("t_refine", r.seconds_refine)
            .num("t_locate", r.seconds_locate)
            .num("refine_clamped", r.n_clamped);
        report_json(j, &rep.simplify, &comp_in, v0);
        certify(j, &coarse);
    } else {
        let mut coarse = CoarseMesh::from_mesh(&mesh, &geometry, &[], CoarseMeshOptions::default())
            .map_err(|e| format!("{e}"))?;
        let rep = simplify(&mut coarse, &so).map_err(|e| format!("{e}"))?;
        report_json(j, &rep, &comp_in, v0);
        certify(j, &coarse);
    }
    Ok(())
}

fn certify(j: &mut Json, coarse: &CoarseMesh) {
    let t = Instant::now();
    match coarse.certify(false) {
        Ok(c) => {
            j.num("certify_ok", 1)
                .num("certify_gauss_bonnet_error", c.gauss_bonnet_error)
                .num("certify_max_angle_error", c.max_angle_error)
                .num("certify_not_delaunay", c.n_not_delaunay);
        }
        Err(e) => {
            j.num("certify_ok", 0).str("certify_error", &e.to_string());
        }
    }
    j.num("t_certify", t.elapsed().as_secs_f64());
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(path) = args.first().cloned() else {
        eprintln!("usage: thingi_simplify MESH.obj [--id ID] [--no-refine] [--signposts] [--target-ratio R]");
        std::process::exit(2);
    };
    let id = args
        .iter()
        .position(|a| a == "--id")
        .and_then(|i| args.get(i + 1))
        .cloned()
        .unwrap_or_default();
    let total = Instant::now();
    let mut j = Json::default();
    j.str("id", &id);
    let outcome = catch_unwind(AssertUnwindSafe(|| run(&path, &args, &mut j)));
    j.num("t_total", total.elapsed().as_secs_f64());
    j.num("peak_rss_mb", peak_rss_kb() / 1024);
    match outcome {
        Ok(Ok(())) => {
            j.str("status", "ok");
        }
        Ok(Err(e)) => {
            j.str("status", "fail").str("error", &e);
        }
        Err(p) => {
            let msg = p
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_default();
            j.str("status", "panic").str("error", &msg);
        }
    }
    println!("{{{}}}", j.0);
}
