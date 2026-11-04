//! The saddle demonstration (`SADDLE-DEMO.md`): our crate's side of the three-way comparison on
//! saddle-dominated analytic shapes whose inputs (points, exact normals, the *exact* Monge
//! coefficients of the true second fundamental form, grids, true SDF, a fine truth mesh) are
//! written once by `runs_tmp/pat-saddle/driver/make_inputs.py` and read by both this example and
//! the driver of the authors' program.
//!
//! ```text
//! saddle_demo verify <root>                 exact coefficients of a sphere and of a torus's equators
//!                                           through our fit (both sign rules): radii, centre, axis, κ
//! saddle_demo case <root> <case> [--threads 8]
//!     (N)     Model::reference() + fit, SignRule::Reference (as the authors) and Geometric
//!     (E-ref) exact coefficients → fit_tori_with(SignRule::Reference) → TorusField::new
//!     (E-geo) the same with SignRule::Geometric
//!     gwn     our winding number (Barill point areas) for the sign, nearest-sample distance
//!     → <root>/cases/<case>/ours/grid64_<m>.npy (f64, 32-NN), grid64_<m>-all.npy (all points),
//!       grid128_<m>.npy (f32), tori_<m>_{centers,axes,major,minor}.npy, diag_<m>.npy (κmin, κmax
//!       of fit_diagnostics, fit convention: convex negative), truth_check.json
//! saddle_demo render <root> <case> <out-dir> [--res 512] [--threads 8]
//!     sphere-traced zero level sets (truth from the fine mesh, N, E-ref, E-geo): a view into the hole,
//!     a cutaway (max(φ, ±y)) and the renderer's cut plane through the axis (φ with isobands)
//! ```
//!
//! The sign rule is always passed explicitly. Every stage runs on an explicit pool (`--threads`).
//! The GPU is never used (CPU device everywhere).

// the shared helpers (`eval_util/util.rs`) are written with `!(x > 0.0)` NaN-catching tests, as in `eval`
#![allow(clippy::neg_cmp_op_on_partial_ord)]

#[path = "eval_util/util.rs"]
mod util;

use ddg_points_as_tori::baselines::{PointAreaOptions, point_areas, winding_numbers, winding_signed_distance};
use ddg_points_as_tori::csg::{Csg, ImplicitField};
use ddg_points_as_tori::fit::fit_tori_with;
use ddg_points_as_tori::render::{Camera, CutPlane, TraceOptions, Tracer};
use ddg_points_as_tori::spatial::{KdTree, KdTreeOptions};
use ddg_points_as_tori::truth::{MeshTruth, SignMode, TruthOptions};
use ddg_points_as_tori::{
    BuildOptions, EigenRule, EvalOptions, EvalOut, FitOptions, Model, Neighborhood, NormalPolicy, PatError, PointCloud, PointsAsTori,
    SignRule, Tori, TorusField, fit_diagnostics,
};
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;
use util::*;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let r = match args.first().map(String::as_str) {
        Some("verify") => verify(&args[1..]),
        Some("case") => case(&args[1..]),
        Some("render") => render(&args[1..]),
        _ => Err("usage: saddle_demo verify|case|render ... (see the file header)".into()),
    };
    if let Err(e) = r {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn opt(a: &[String], key: &str, default: usize) -> R<usize> {
    match a.iter().position(|s| s == key) {
        Some(i) => a.get(i + 1).ok_or(format!("{key} needs a value"))?.parse().map_err(|_| format!("{key}: not a number")),
        None => Ok(default),
    }
}

fn read_coeffs(path: &Path) -> R<Vec<[f32; 6]>> {
    let (s, d) = read_flat(path)?;
    if s.len() != 2 || s[1] != 6 {
        return Err(format!("{}: expected (n, 6), got {s:?}", path.display()));
    }
    // f32 arrays are returned bit-exact as f64, so the cast back is exact
    Ok(d.chunks_exact(6).map(|c| std::array::from_fn(|i| c[i] as f32)).collect())
}

fn load_cloud(dir: &Path) -> R<PointCloud> {
    let p = read3(&dir.join("points.npy"))?;
    let n = read3(&dir.join("normals.npy"))?;
    Ok(PointCloud::new(p, n, NormalPolicy::RequireUnit { tol: 1e-9 }).map_err(err)?.0)
}

fn fit_opts(rule: SignRule, threads: usize) -> FitOptions {
    // explicit: sign rule given, eigenvector rule as trained and as the authors' program
    FitOptions { sign_rule: rule, eigen_rule: EigenRule::Reference, threads: Some(threads), ..Default::default() }
}

fn eval_opts(all: bool, threads: usize) -> EvalOptions {
    let neighborhood = if all { Neighborhood::All } else { Neighborhood::Knn { k: 32 } };
    EvalOptions { neighborhood, threads: Some(threads), ..Default::default() }
}

/// The tori of `coeffs` fitted with `rule` (no network): the exact-coefficient field.
fn exact_field(cloud: &PointCloud, coeffs: &[[f32; 6]], rule: SignRule, all: bool, threads: usize) -> R<TorusField> {
    let sigma = vec![f64::NAN; cloud.len()];
    let (tori, _) = fit_tori_with(cloud, coeffs, &sigma, 0, &fit_opts(rule, threads)).map_err(err)?;
    Ok(TorusField::new(cloud, tori, eval_opts(all, threads)).map_err(err)?.with_sign_rule(Some(rule)))
}

/// The network's field with `rule` (coefficients predicted once, refitted per rule).
fn network_fields(cloud: &PointCloud, threads: usize) -> R<(Tori, TorusField, TorusField)> {
    let model = Model::reference().map_err(err)?;
    let mut bo = BuildOptions::default().with_threads(Some(threads));
    bo.fit = fit_opts(SignRule::Reference, threads);
    bo.eval = eval_opts(false, threads);
    let (fref, _) = PointsAsTori::build(cloud, &model, &bo).map_err(err)?;
    let t = fref.tori().clone();
    let (tg, _) = fit_tori_with(cloud, &t.coefficients, &t.sigma, 0, &fit_opts(SignRule::Geometric, threads)).map_err(err)?;
    let fgeo = TorusField::new(cloud, tg, eval_opts(false, threads)).map_err(err)?.with_sign_rule(Some(SignRule::Geometric));
    Ok((t, fref, fgeo))
}

fn grid_of(axis: &[f64]) -> Vec<[f64; 3]> {
    let mut g = Vec::with_capacity(axis.len().pow(3));
    for &x in axis {
        for &y in axis {
            for &z in axis {
                g.push([x, y, z]);
            }
        }
    }
    g
}

fn eval_grid(f: &TorusField, xs: &[[f64; 3]]) -> R<Vec<f64>> {
    let mut phi = vec![0.0; xs.len()];
    f.evaluate(xs, EvalOut { distance: Some(&mut phi), ..Default::default() }).map_err(err)?;
    Ok(phi)
}

fn save_tori(dir: &Path, tag: &str, t: &Tori) -> R<()> {
    let c: Vec<[f64; 3]> = t.tori.iter().map(|t| t.center).collect();
    let a: Vec<[f64; 3]> = t.tori.iter().map(|t| t.axis).collect();
    write3(&dir.join(format!("tori_{tag}_centers.npy")), &c)?;
    write3(&dir.join(format!("tori_{tag}_axes.npy")), &a)?;
    let n = t.tori.len();
    write_npy_f64(&dir.join(format!("tori_{tag}_major.npy")), &[n], &t.tori.iter().map(|t| t.major).collect::<Vec<_>>())?;
    write_npy_f64(&dir.join(format!("tori_{tag}_minor.npy")), &[n], &t.tori.iter().map(|t| t.minor).collect::<Vec<_>>())
}

fn save_diag(dir: &Path, tag: &str, cloud: &PointCloud, coeffs: &[[f32; 6]], rule: SignRule) -> R<()> {
    let o = fit_opts(rule, 1);
    let mut v = Vec::with_capacity(2 * cloud.len());
    for i in 0..cloud.len() {
        let d = fit_diagnostics(cloud, coeffs, i, &o).map_err(err)?;
        v.push(d.kappa_min);
        v.push(d.kappa_max);
    }
    write_npy_f64(&dir.join(format!("diag_{tag}.npy")), &[cloud.len(), 2], &v)
}

// ------------------------------------------------------------------ verify

fn verify(a: &[String]) -> R<()> {
    let root = PathBuf::from(a.first().ok_or("verify <root>")?);
    let d = root.join("verify");
    let cloud = load_cloud(&d)?;
    let coeffs = read_coeffs(&d.join("coeffs_exact.npy"))?;
    let (_, kt) = read_flat(&d.join("kappa_true.npy"))?;
    let expected = std::fs::read_to_string(d.join("expected.json")).map_err(err)?;
    // the expected values are written by make_inputs.py in a fixed order: 8 sphere points (centre 0,
    // R 0, r 0.5), 6 rotated outer-equator points (centre 0, R .54, r .36), 6 rotated inner-equator
    // points (R .54, r .18), 3 axis-aligned outer-equator points (degenerate pair, reported only)
    if !expected.contains("torus-inner-equator") {
        return Err("verify/expected.json: unexpected content".into());
    }
    let mut j = Json::new();
    let mut lines = Vec::new();
    for (rule, eig, tag) in [
        (SignRule::Reference, EigenRule::Reference, "reference"),
        (SignRule::Geometric, EigenRule::Reference, "geometric"),
        (SignRule::Geometric, EigenRule::Robust, "geometric-robust"),
    ] {
        let o = FitOptions { sign_rule: rule, eigen_rule: eig, threads: Some(1), ..Default::default() };
        let sigma = vec![f64::NAN; cloud.len()];
        let (tori, _) = fit_tori_with(&cloud, &coeffs, &sigma, 0, &o).map_err(err)?;
        let mut rows = Vec::new();
        for (i, t) in tori.tori.iter().enumerate() {
            let dg = fit_diagnostics(&cloud, &coeffs, i, &o).map_err(err)?;
            // fit convention is convex negative; the truth convex positive
            let ek = ((-dg.kappa_min - kt[2 * i]).abs()).max((-dg.kappa_max - kt[2 * i + 1]).abs());
            let mut r = Json::new();
            r.int("i", i as i64)
                .raw("center", &format!("[{:.9},{:.9},{:.9}]", t.center[0], t.center[1], t.center[2]))
                .raw("axis", &format!("[{:.9},{:.9},{:.9}]", t.axis[0], t.axis[1], t.axis[2]))
                .num("major", t.major)
                .num("minor", t.minor)
                .num("kappa_err", ek);
            rows.push(r.finish());
            let p = cloud.positions()[i];
            let n = cloud.normals()[i];
            // osculation: g(p) = 0 and ∇g(p) = n (the torus is on the correct side)
            let (g, gr, _) = ddg_points_as_tori::sdf::torus_sdf_gradient(t, p);
            let dn = gr[0] * n[0] + gr[1] * n[1] + gr[2] * n[2];
            lines.push(format!(
                "{tag:17} {i:2}: R {:.6} r {:+.6} centre ({:+.6},{:+.6},{:+.6}) |κ err| {:.1e} g(p) {:+.1e} ∇g·n {:+.3}",
                t.major, t.minor, t.center[0], t.center[1], t.center[2], ek, g, dn
            ));
        }
        j.raw(tag, &jarr(&rows));
    }
    std::fs::write(d.join("verify_ours.json"), j.finish()).map_err(err)?;
    for l in &lines {
        println!("{l}");
    }
    Ok(())
}

// ------------------------------------------------------------------ case

fn case(a: &[String]) -> R<()> {
    let root = PathBuf::from(a.first().ok_or("case <root> <case>")?);
    let name = a.get(1).ok_or("case <root> <case>")?;
    let th = opt(a, "--threads", 8)?;
    let dir = root.join("cases").join(name);
    let out = dir.join("ours");
    std::fs::create_dir_all(&out).map_err(err)?;
    let t0 = Instant::now();
    let cloud = load_cloud(&dir)?;
    let exact = read_coeffs(&dir.join("coeffs_exact.npy"))?;
    let (_, ax64) = read_flat(&root.join("grid_axis64.npy"))?;
    let (_, ax128) = read_flat(&root.join("grid_axis128.npy"))?;
    let (g64, g128) = (grid_of(&ax64), grid_of(&ax128));
    // the plane y = 0 on 384 x 384 points of [-1, 1]^2, index [ix, iz] (as the authors' driver)
    let sl: Vec<f64> = (0..384).map(|i| -1.0 + 2.0 * i as f64 / 383.0).collect();
    let gs: Vec<[f64; 3]> = sl.iter().flat_map(|&x| sl.iter().map(move |&z| [x, 0.0, z])).collect();
    let n = cloud.len();
    let mut log = Json::new();
    log.str("case", name).int("n", n as i64).int("threads", th as i64);

    let save64 = |tag: &str, v: &[f64]| write_npy_f64(&out.join(format!("grid64_{tag}.npy")), &[64, 64, 64], v);
    let save128 = |tag: &str, v: &[f64]| write_npy_f32(&out.join(format!("grid128_{tag}.npy")), &[128, 128, 128], v);
    let save_sl = |tag: &str, v: &[f64]| write_npy_f64(&out.join(format!("slice_{tag}.npy")), &[384, 384], v);

    // (N) with both rules
    let (tn, fnref, fngeo) = network_fields(&cloud, th)?;
    save_tori(&out, "Nref", fnref.tori())?;
    save_tori(&out, "Ngeo", fngeo.tori())?;
    write_npy_f64(&out.join("net_sigma.npy"), &[n], &tn.sigma)?;
    let flat: Vec<f64> = tn.coefficients.iter().flatten().map(|&x| x as f64).collect();
    write_npy_f32(&out.join("net_coeffs.npy"), &[n, 6], &flat)?;
    save_diag(&out, "N", &cloud, &tn.coefficients, SignRule::Reference)?;
    save64("Nref", &eval_grid(&fnref, &g64)?)?;
    save64("Ngeo", &eval_grid(&fngeo, &g64)?)?;
    save128("Nref", &eval_grid(&fnref, &g128)?)?;
    save_sl("Nref", &eval_grid(&fnref, &gs)?)?;
    save_sl("Ngeo", &eval_grid(&fngeo, &gs)?)?;
    let all = TorusField::new(&cloud, fnref.tori().clone(), eval_opts(true, th)).map_err(err)?;
    save64("Nref-all", &eval_grid(&all, &g64)?)?;

    // (E-ref), (E-geo)
    save_diag(&out, "E", &cloud, &exact, SignRule::Reference)?;
    for (rule, tag) in [(SignRule::Reference, "Eref"), (SignRule::Geometric, "Egeo")] {
        let f = exact_field(&cloud, &exact, rule, false, th)?;
        save_tori(&out, tag, f.tori())?;
        save64(tag, &eval_grid(&f, &g64)?)?;
        save128(tag, &eval_grid(&f, &g128)?)?;
        save_sl(tag, &eval_grid(&f, &gs)?)?;
        let fa = exact_field(&cloud, &exact, rule, true, th)?;
        save64(&format!("{tag}-all"), &eval_grid(&fa, &g64)?)?;
    }

    // GWN sign + nearest-sample distance
    let (areas, _) = point_areas(&cloud, &PointAreaOptions { threads: Some(th), ..Default::default() }).map_err(err)?;
    let tree = KdTree::new(cloud.positions(), KdTreeOptions::default()).map_err(err)?;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(th).build().map_err(err)?;
    for (g, n_) in [(&g64, 64usize), (&g128, 128), (&gs, 384)] {
        let mut w = vec![0.0; g.len()];
        winding_numbers(&cloud, &areas, g, &mut w, Some(th)).map_err(err)?;
        let near: Vec<f64> = pool.install(|| {
            g.par_chunks(1024)
                .flat_map_iter(|c| {
                    let mut buf = Vec::new();
                    c.iter()
                        .map(|x| {
                            tree.knn(*x, 1, &mut buf);
                            buf[0].d2.sqrt()
                        })
                        .collect::<Vec<_>>()
                })
                .collect()
        });
        let s: Vec<f64> = w.iter().zip(&near).map(|(&w, &d)| winding_signed_distance(d, w, 0.5)).collect();
        match n_ {
            64 => save64("gwn", &s)?,
            128 => save128("gwn", &s)?,
            _ => save_sl("gwn", &s)?,
        }
    }
    log.num("area_total", areas.iter().sum());

    // cross-check of the fine truth mesh (used by the renders) against the analytic truth grid
    let v = read3(&dir.join("mesh_vertices.npy"))?;
    let f = read_faces(&dir.join("mesh_faces.npy"))?;
    let (mt, rep) = MeshTruth::new(&v, &f, &TruthOptions { sign: SignMode::Exact, threads: Some(th), ..Default::default() }).map_err(err)?;
    let (_, t64) = read_flat(&dir.join("truth_sdf64.npy"))?;
    // every 7th grid point (37k queries; the exact winding number is a sum over all faces)
    let sub: Vec<usize> = (0..g64.len()).step_by(7).collect();
    let xs: Vec<[f64; 3]> = sub.iter().map(|&i| g64[i]).collect();
    let mut s = vec![0.0; xs.len()];
    mt.evaluate(&xs, &mut s, None).map_err(err)?;
    let dmax = sub.iter().zip(&s).map(|(&i, v)| (v - t64[i]).abs()).fold(0.0, f64::max);
    let sign_bad = sub.iter().zip(&s).filter(|(i, v)| (**v < 0.0) != (t64[**i] < 0.0) && t64[**i].abs() > 1e-3).count();
    log.int("mesh_faces", rep.faces as i64).num("mesh_volume", rep.signed_volume).num("mesh_vs_analytic_max", dmax).int("mesh_vs_analytic_sign", sign_bad as i64);
    log.num("seconds", t0.elapsed().as_secs_f64());
    let s = log.finish();
    std::fs::write(out.join("log.json"), &s).map_err(err)?;
    println!("{s}");
    Ok(())
}

// ------------------------------------------------------------------ render

/// A field with a fixed bounding box (so every method is framed identically).
struct Framed {
    f: Arc<dyn ImplicitField + Send + Sync>,
    bb: [[f64; 3]; 2],
}

impl ImplicitField for Framed {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        self.f.distance(x)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        self.f.distance_gradient(x)
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        self.bb
    }
    fn distances(&self, xs: &[[f64; 3]], out: &mut [f64]) -> Result<(), PatError> {
        self.f.distances(xs, out)
    }
}

/// The true surface from the fine mesh: exact signed distance, gradient `±(x − closest)/d`.
struct MeshField(MeshTruth);

impl ImplicitField for MeshField {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        self.0.signed_distance(x)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let s = self.0.signed_distance(x)?;
        let (d, c, _) = self.0.closest(x);
        if !(d > 1e-12) {
            return Ok((s, [0.0, 0.0, 1.0]));
        }
        let k = s.signum() / d;
        Ok((s, [(x[0] - c[0]) * k, (x[1] - c[1]) * k, (x[2] - c[2]) * k]))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        [[-1.0; 3], [1.0; 3]]
    }
}

/// `φ = s·(x·e_axis) + c` (a half-space; with `max` it cuts the solid away where `s·x_axis > −c`).
struct HalfSpace {
    axis: usize,
    s: f64,
    c: f64,
}

impl ImplicitField for HalfSpace {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        Ok(self.s * x[self.axis] + self.c)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let mut g = [0.0; 3];
        g[self.axis] = self.s;
        Ok((self.s * x[self.axis] + self.c, g))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        [[-1.0; 3], [1.0; 3]]
    }
}

/// The field with everything on the camera side of the plane `y = 0` removed (a constant positive
/// value for `y > 1e-9`), so that the renderer's cut plane at `y = 0` shows the field's own `φ`
/// everywhere on the plane (the tracer stops at the plane, which is then drawn).
struct FrontCut(Arc<dyn ImplicitField + Send + Sync>);

impl ImplicitField for FrontCut {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        if x[1] > 1e-9 { Ok(1.0) } else { self.0.distance(x) }
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        if x[1] > 1e-9 { Ok((1.0, [0.0, 1.0, 0.0])) } else { self.0.distance_gradient(x) }
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        self.0.bbox()
    }
    fn distances(&self, xs: &[[f64; 3]], out: &mut [f64]) -> Result<(), PatError> {
        self.0.distances(xs, out)?;
        for (x, o) in xs.iter().zip(out.iter_mut()) {
            if x[1] > 1e-9 {
                *o = 1.0;
            }
        }
        Ok(())
    }
}

fn render(a: &[String]) -> R<()> {
    let root = PathBuf::from(a.first().ok_or("render <root> <case> <out-dir>")?);
    let name = a.get(1).ok_or("render <root> <case> <out-dir>")?;
    let outdir = PathBuf::from(a.get(2).ok_or("render <root> <case> <out-dir>")?);
    let th = opt(a, "--threads", 8)?;
    let res = opt(a, "--res", 512)?;
    std::fs::create_dir_all(&outdir).map_err(err)?;
    let dir = root.join("cases").join(name);
    let cloud = load_cloud(&dir)?;
    let exact = read_coeffs(&dir.join("coeffs_exact.npy"))?;
    let v = read3(&dir.join("mesh_vertices.npy"))?;
    let f = read_faces(&dir.join("mesh_faces.npy"))?;
    let (mt, _) = MeshTruth::new(&v, &f, &TruthOptions { sign: SignMode::Fast { beta: 2.0 }, threads: Some(th), ..Default::default() }).map_err(err)?;
    let (_, fnref, _) = network_fields(&cloud, th)?;
    let methods: Vec<(&str, Arc<dyn ImplicitField + Send + Sync>)> = vec![
        ("truth", Arc::new(MeshField(mt))),
        ("N", Arc::new(fnref)),
        ("Eref", Arc::new(exact_field(&cloud, &exact, SignRule::Reference, false, th)?)),
        ("Egeo", Arc::new(exact_field(&cloud, &exact, SignRule::Geometric, false, th)?)),
    ];
    let bb = [[-1.0; 3], [1.0; 3]];
    let spool = name.starts_with("spool");
    for (tag, field) in &methods {
        let framed = |f: Arc<dyn ImplicitField + Send + Sync>| Framed { f, bb };
        // 1. a view of the whole surface (into the hole for the tori, side view for the spool)
        // 2. a cutaway: the half y > 0 (camera side) removed, seen from +y above: the tube cross-sections
        // 3. the renderer's cut plane y = 0 coloured by the method's φ with isobands every 0.05, seen
        //    from +y with the camera side removed (FrontCut)
        let views: [(&str, f64, f64, Option<HalfSpace>, bool); 3] = [
            ("view", if spool { 30.0 } else { 20.0 }, if spool { 15.0 } else { 55.0 }, None, false),
            ("cutaway", 90.0, 25.0, Some(HalfSpace { axis: 1, s: 1.0, c: 0.0 }), false),
            ("cutplane", 90.0, 0.0, None, true),
        ];
        for (vname, az, el, cut, plane) in views {
            let scene: Box<dyn ImplicitField> = match (cut, plane) {
                (Some(h), _) => Box::new(framed(Arc::new(Csg::Leaf(field.clone()).intersect(Csg::leaf(h))))),
                (None, true) => Box::new(framed(Arc::new(FrontCut(field.clone())))),
                (None, false) => Box::new(framed(field.clone())),
            };
            let mut o = TraceOptions { threads: Some(th), ..Default::default() };
            if let Some(g) = o.ground.as_mut() {
                g.up = [0.0, 0.0, 1.0];
            }
            if vname != "view" {
                o.ground = None;
            }
            if plane {
                o.cut_plane = Some(CutPlane { axis: 1, offset: 0.0, range: Some(0.5), isoband_frequency: Some(2.0 * std::f64::consts::PI / 0.05), ..Default::default() });
            }
            let cam = Camera::orbit([0.0; 3], 12.0, az.to_radians(), el.to_radians(), 2, 11.0, res, res).map_err(err)?;
            let tracer = Tracer::new(o).map_err(err)?;
            let t = Instant::now();
            let (img, rep) = tracer.render(scene.as_ref(), &cam).map_err(err)?;
            let path = outdir.join(format!("{name}_{tag}_{vname}.png"));
            img.write_png(&path, false).map_err(err)?;
            println!(
                "{}: {:.1}s hits {} misses {} step-limit {} non-finite {}",
                path.display(),
                t.elapsed().as_secs_f64(),
                rep.hits,
                rep.misses,
                rep.step_limit,
                rep.non_finite
            );
        }
    }
    Ok(())
}
