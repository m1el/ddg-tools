//! WP10, SPEC §10.2: the Thingi10k evaluation of one mesh with one representation.
//!
//! ```text
//! thingi_eval <mesh.obj> --backend integer|sp-ruppert|sp-chew [--id ID]
//! ```
//!
//! One process, one mesh, one backend; prints one JSON line on stdout (also on failure). The
//! driver `scripts/ic-eval.sh` runs it in parallel with a time cap and a memory limit.
//!
//! Mesh: edge lengths from the positions, intrinsic mollification with `ε = 1e-5 h̄` (as the
//! paper, Sec. 3.9), then
//!
//! * **IDT stage**: flips to Delaunay; the integers are certified (SPEC §2.4), the signposts'
//!   edges are traced over the input (every trace passes its endpoint check); the common
//!   subdivision is built (integers: `common_subdivision`; signposts: `overlay` of the traces,
//!   IT-SPEC N12) and checked: fan-triangulated it must be a valid `TriMesh` with the Euler
//!   characteristic, boundary edges and components of the input.
//! * **Refinement stage**: refinement to 25° from the IDT (integers: Chew's rules, the paper's;
//!   signposts: Ruppert's rules, the signposts' default, or Chew's), certificate / traces, the
//!   same subdivision and the same checks.
//!
//! Success of a stage is "the common subdivision's connectivity is recovered" (Table 1). The
//! signposts are run on the mesh as it is, with its boundary if any (the table of the
//! results separates closed meshes, which is what the paper evaluated).

use std::fmt::Write as _;
use std::time::Instant;

use ddg_core::io::parse_obj;
use ddg_core::mollify::mollify_edge_lengths;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::{
    InitOptions, IntegerRefineOptions, IntegerTriangulation, PAPER_MOLLIFY_EPSILON,
};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, IntrinsicTriangulation, Overlay, RefineOptions,
    RefineReport, RefineRules, Triangulation,
};

/// A JSON object under construction.
#[derive(Default)]
struct Json(String);

impl Json {
    fn sep(&mut self) {
        if !self.0.is_empty() {
            self.0.push(',');
        }
    }
    fn num(&mut self, k: &str, v: impl std::fmt::Display) -> &mut Self {
        self.sep();
        let s = v.to_string();
        let s = if s == "inf" || s == "-inf" || s == "NaN" {
            "null".to_string()
        } else {
            s
        };
        let _ = write!(self.0, "\"{k}\":{s}");
        self
    }
    fn str(&mut self, k: &str, v: &str) -> &mut Self {
        self.sep();
        let esc: String = v
            .chars()
            .flat_map(|c| match c {
                '"' => vec!['\\', '"'],
                '\\' => vec!['\\', '\\'],
                '\n' | '\t' | '\r' => vec![' '],
                c => vec![c],
            })
            .collect();
        let esc: String = esc.chars().take(400).collect();
        let _ = write!(self.0, "\"{k}\":\"{esc}\"");
        self
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

/// The N12-style checks of a subdivision: a valid fan-triangulated `TriMesh` with the topology
/// of the input. Returns `(V, E, F)` of the polygon mesh, or the reason.
fn check_subdivision(
    sub: &Overlay,
    input: &TriMesh,
) -> Result<(usize, usize, usize, usize), String> {
    let nv = sub.vertices.len();
    let nf = sub.n_faces();
    for f in 0..nf {
        let c = sub.face(f);
        if !(3..=6).contains(&c.len()) {
            return Err(format!("piece {f} has {} corners", c.len()));
        }
    }
    // Edges of the polygon mesh: directed keys, sorted.
    let mut keys: Vec<u64> = Vec::with_capacity(sub.corners.len());
    for f in 0..nf {
        let c = sub.face(f);
        for k in 0..c.len() {
            let (a, b) = (c[k], c[(k + 1) % c.len()]);
            keys.push((u64::from(a) << 32) | u64::from(b));
        }
    }
    keys.sort_unstable();
    if let Some(w) = keys.windows(2).find(|w| w[0] == w[1]) {
        return Err(format!(
            "directed edge ({}, {}) occurs twice",
            w[0] >> 32,
            w[0] & 0xffff_ffff
        ));
    }
    let ne = {
        let mut und: Vec<u64> = keys
            .iter()
            .map(|&k| {
                let (a, b) = (k >> 32, k & 0xffff_ffff);
                (a.min(b) << 32) | a.max(b)
            })
            .collect();
        und.sort_unstable();
        und.dedup();
        und.len()
    };
    drop(keys);
    let (tris, _) = sub.triangles();
    let mesh = TriMesh::new(nv, tris).map_err(|e| format!("fan triangulation: {e}"))?;
    if mesh.euler_characteristic() != input.euler_characteristic() {
        return Err(format!(
            "Euler characteristic {} instead of {}",
            mesh.euler_characteristic(),
            input.euler_characteristic()
        ));
    }
    if mesh.n_components() != input.n_components() {
        return Err(format!(
            "{} components instead of {}",
            mesh.n_components(),
            input.n_components()
        ));
    }
    if mesh.boundary_loops().len() != input.boundary_loops().len() {
        return Err(format!(
            "{} boundary loops instead of {}",
            mesh.boundary_loops().len(),
            input.boundary_loops().len()
        ));
    }
    Ok((nv, ne, nf, sub.corners.len() - 2 * nf))
}

#[derive(Default)]
struct Stage {
    ok: bool,
    error: String,
    seconds: f64,
}

fn record_subdivision(
    j: &mut Json,
    prefix: &str,
    r: Result<(usize, usize, usize, usize), String>,
) -> Result<(), String> {
    match r {
        Ok((v, e, f, tf)) => {
            j.num(&format!("{prefix}_sub_v"), v)
                .num(&format!("{prefix}_sub_e"), e)
                .num(&format!("{prefix}_sub_f"), f)
                .num(&format!("{prefix}_sub_tf"), tf);
            Ok(())
        }
        Err(m) => Err(m),
    }
}

fn refine_fields(j: &mut Json, r: &RefineReport) {
    j.num("inserted", r.n_inserted)
        .num("removals", r.n_removals)
        .num("removals_refused", r.n_removals_refused)
        .num("trace_failures", r.n_trace_failures)
        .num("min_angle_deg", r.min_angle.to_degrees())
        .num("min_angle_valid_deg", r.min_angle_valid.to_degrees())
        .num(
            "min_angle_reachable_deg",
            r.min_angle_reachable.to_degrees(),
        )
        .num("n_below_bound", r.below_bound.len())
        .num("n_not_delaunay", r.n_not_delaunay)
        .num("n_encroached", r.n_encroached)
        .num("refine_flips", r.n_flips)
        .num(
            "insertion_limit",
            u8::from(matches!(
                r.status,
                ddg_intrinsic::RefineStatus::InsertionLimit
            )),
        );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let path = args.first().expect("mesh path").clone();
    let mut backend = String::from("integer");
    let mut id = String::new();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--backend" => backend = args[i + 1].clone(),
            "--id" => id = args[i + 1].clone(),
            other => panic!("unknown argument {other}"),
        }
        i += 2;
    }
    let total = Instant::now();
    let mut j = Json::default();
    j.str("id", &id).str("backend", &backend);
    let outcome = run(&path, &backend, &mut j);
    j.num("seconds_total", total.elapsed().as_secs_f64());
    j.num("peak_rss_mb", peak_rss_kb() / 1024);
    match outcome {
        Ok(()) => {
            j.str("status", "ok");
        }
        Err(e) => {
            j.str("status", "fail").str("error", &e);
        }
    }
    println!("{{{}}}", j.0);
}

fn run(path: &str, backend: &str, j: &mut Json) -> Result<(), String> {
    let t = Instant::now();
    let text = std::fs::read(path).map_err(|e| format!("read: {e}"))?;
    let obj = parse_obj(&String::from_utf8_lossy(&text)).map_err(|e| format!("parse: {e}"))?;
    let faces = obj.triangles().ok_or("not all triangles")?;
    let mesh = TriMesh::new(obj.positions.len(), faces).map_err(|e| format!("mesh: {e}"))?;
    let p = &obj.positions;
    let lengths: Vec<f64> = mesh
        .edges()
        .iter()
        .map(|&[a, b]| {
            let (a, b) = (p[a as usize], p[b as usize]);
            ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
        })
        .collect();
    drop(text);
    j.num("v0", mesh.n_vertices())
        .num("f0", mesh.n_faces())
        .num("boundary_edges", mesh.n_boundary_edges())
        .num("components", mesh.n_components())
        .num("chi", mesh.euler_characteristic())
        .num("seconds_load", t.elapsed().as_secs_f64());
    let delaunay = DelaunayOptions::default();
    // The paper's benchmark runs with --refineMaxInsertions=0: no limit.
    let refine_options = RefineOptions {
        max_insertions: ddg_intrinsic::InsertionLimit::Unlimited,
        ..RefineOptions::default()
    };
    let mut idt = Stage::default();
    let mut refn = Stage::default();

    if backend == "integer" {
        let t = Instant::now();
        let (mut tri, init) = IntegerTriangulation::with_options(
            &mesh,
            &lengths,
            &InitOptions {
                mollify: Some(PAPER_MOLLIFY_EPSILON),
            },
        )
        .map_err(|e| format!("init: {e}"))?;
        if let Some(m) = &init.mollified {
            j.num("mollify_epsilon", m.epsilon)
                .num("mollify_faces_below", m.n_faces_below_margin);
        }
        j.num("seconds_init", t.elapsed().as_secs_f64());
        // IDT.
        let t = Instant::now();
        let d = tri
            .flip_to_delaunay_certified(&delaunay)
            .map_err(|e| format!("idt: {e}"))?;
        idt.seconds = t.elapsed().as_secs_f64();
        dump_idt(tri.edge_lengths());
        j.num("idt_seconds", idt.seconds)
            .num("idt_flips", d.0.n_flips)
            .num("idt_not_delaunay", d.0.n_not_delaunay)
            .num("idt_certificate_ok", 1)
            .num("idt_sum_nplus", tri.total_crossings());
        let t = Instant::now();
        match tri.common_subdivision() {
            Ok(sub) => {
                j.num("idt_sub_seconds", t.elapsed().as_secs_f64())
                    .num("idt_clamps", sub.report.strip.clamped_crossings);
                match record_subdivision(
                    j,
                    "idt",
                    check_subdivision(&sub.overlay, tri.input().mesh()),
                ) {
                    Ok(()) => idt.ok = true,
                    Err(m) => idt.error = m,
                }
            }
            Err(e) => idt.error = format!("subdivision: {e}"),
        }
        j.num("idt_ok", u8::from(idt.ok))
            .str("idt_error", &idt.error);
        // Refinement.
        let t = Instant::now();
        let options = IntegerRefineOptions {
            refine: refine_options,
            rules: RefineRules::Chew,
            ..IntegerRefineOptions::default()
        };
        match tri.refine(&options) {
            Ok(r) => {
                j.num("refine_seconds", t.elapsed().as_secs_f64())
                    .num("refine_certificate_ok", u8::from(r.certificate.is_some()))
                    .num("clamped_crossings", r.clamped_crossings)
                    .num("face_split_refused", r.n_recover_failed)
                    .num("near_arc", r.n_near_arc)
                    .num("places_moved", r.n_places_moved);
                refine_fields(j, &r.refine);
                j.num("v1", tri.complex().n_vertices())
                    .num("refine_sum_nplus", tri.total_crossings());
                let t = Instant::now();
                match tri.common_subdivision() {
                    Ok(sub) => {
                        j.num("refine_sub_seconds", t.elapsed().as_secs_f64())
                            .num("refine_sub_clamps", sub.report.strip.clamped_crossings);
                        match record_subdivision(
                            j,
                            "refine",
                            check_subdivision(&sub.overlay, tri.input().mesh()),
                        ) {
                            Ok(()) => refn.ok = true,
                            Err(m) => refn.error = m,
                        }
                    }
                    Err(e) => refn.error = format!("subdivision: {e}"),
                }
            }
            Err(e) => refn.error = format!("refine: {e}"),
        }
    } else {
        let rules = match backend {
            "sp-ruppert" => RefineRules::Ruppert,
            "sp-chew" => RefineRules::Chew,
            other => return Err(format!("unknown backend {other}")),
        };
        let t = Instant::now();
        let m = mollify_edge_lengths(&mesh, &lengths, PAPER_MOLLIFY_EPSILON)
            .map_err(|e| format!("mollify: {e}"))?;
        j.num("mollify_epsilon", m.epsilon)
            .num("mollify_faces_below", m.n_faces_below_margin);
        let geometry = Geometry::from_edge_lengths(&mesh, m.edge_lengths)
            .map_err(|e| format!("geometry: {e}"))?;
        let mut tri =
            IntrinsicTriangulation::new(&mesh, &geometry).map_err(|e| format!("init: {e}"))?;
        j.num("seconds_init", t.elapsed().as_secs_f64());
        let trace_options = EdgeTraceOptions::default();
        let t = Instant::now();
        let d = tri
            .flip_to_delaunay(&delaunay)
            .map_err(|e| format!("idt: {e}"))?;
        idt.seconds = t.elapsed().as_secs_f64();
        dump_idt(tri.edge_lengths());
        j.num("idt_seconds", idt.seconds)
            .num("idt_flips", d.n_flips)
            .num("idt_not_delaunay", d.n_not_delaunay);
        let t = Instant::now();
        let stage = |tri: &IntrinsicTriangulation| -> Result<(Overlay, usize), String> {
            let traces = tri
                .trace_edges(&trace_options)
                .map_err(|e| format!("traces: {e}"))?;
            let n: usize = (0..traces.n_edges())
                .map(|e| traces.crossings(e).len())
                .sum();
            let o = tri.overlay(&traces).map_err(|e| format!("overlay: {e}"))?;
            Ok((o, n))
        };
        match stage(&tri) {
            Ok((o, n)) => {
                j.num("idt_sub_seconds", t.elapsed().as_secs_f64())
                    .num("idt_sum_nplus", n);
                match record_subdivision(j, "idt", check_subdivision(&o, tri.input().mesh())) {
                    Ok(()) => idt.ok = true,
                    Err(m) => idt.error = m,
                }
            }
            Err(e) => idt.error = e,
        }
        j.num("idt_ok", u8::from(idt.ok))
            .str("idt_error", &idt.error);
        let t = Instant::now();
        match tri.refine_with_rules(&refine_options, rules) {
            Ok(r) => {
                j.num("refine_seconds", t.elapsed().as_secs_f64())
                    .num("snapped", r.n_snapped);
                refine_fields(j, &r);
                j.num("v1", tri.n_vertices());
                let t = Instant::now();
                match stage(&tri) {
                    Ok((o, n)) => {
                        j.num("refine_sub_seconds", t.elapsed().as_secs_f64())
                            .num("refine_sum_nplus", n);
                        match record_subdivision(
                            j,
                            "refine",
                            check_subdivision(&o, tri.input().mesh()),
                        ) {
                            Ok(()) => refn.ok = true,
                            Err(m) => refn.error = m,
                        }
                    }
                    Err(e) => refn.error = e,
                }
            }
            Err(e) => refn.error = format!("refine: {e}"),
        }
    }
    j.num("refine_ok", u8::from(refn.ok))
        .str("refine_error", &refn.error);
    Ok(())
}

/// With `IC_DUMP_IDT=path` in the environment: the sorted edge lengths of the intrinsic
/// Delaunay triangulation, one per line (compared with the reference demo's `faceLengths.dmat`
/// as sets, SPEC 10.3).
fn dump_idt(lengths: &[f64]) {
    if let Ok(path) = std::env::var("IC_DUMP_IDT") {
        let mut l = lengths.to_vec();
        l.sort_by(f64::total_cmp);
        let mut s = String::with_capacity(l.len() * 24);
        for x in l {
            let _ = writeln!(s, "{x:e}");
        }
        let _ = std::fs::write(path, s);
    }
}
