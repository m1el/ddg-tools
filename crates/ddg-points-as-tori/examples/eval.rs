//! WP14 evaluation harness (SPEC.md §16.2–16.4, §18; EVALUATION.md, COMPARISON.md).
//!
//! ```text
//! eval prepare <dataset> <root> [--points N] [--threads T]
//!     datasets: abc blobs thingi analytic saddle analyticsad ref reffull orient
//!     writes <root>/<dataset>[<N>]/<case>/{points,normals}.npy, meta.txt and the truth mesh
//! eval case <case-dir> [--threads T] [--sweep] [--retrained label=path[:geometric]]...
//!     64³ grid errors of every method; writes <case-dir>/result.json and slices.npy
//! eval continuity <case-dir> [--threads T] [--res 64,128,256] [--top 40]
//!     jump statistic, error, eikonality and zero-level-set components per evaluation rule;
//!     writes <case-dir>/continuity.json
//! eval orient <case-dir> [--threads T] [--flip 0.3] [--seed S]
//!     Fig. 27 protocol at 10 and 30 iterations; writes <case-dir>/orient.json
//! ```
//!
//! Every method runs on the CPU (the GPU belongs to the training run). Bad input is reported as an
//! error message and a non-zero exit code, never a panic.

#[path = "eval_util/util.rs"]
mod util;

use ddg_points_as_tori::analytic::{Instance, PAPER_KINDS, Solid, saddle_set};
use ddg_points_as_tori::apps::orient::{OrientOptions, orient_normals};
use ddg_points_as_tori::baselines::{
    KernelField, KernelMethod, KernelOptions, PointAreaOptions, point_areas, regularized_winding_numbers, winding_numbers,
    winding_signed_distance,
};
use ddg_points_as_tori::blend::BlendWant;
use ddg_points_as_tori::csg::ImplicitField;
use ddg_points_as_tori::generate::blob::{BlobOptions, generate_blob};
use ddg_points_as_tori::generate::rng::Rng;
use ddg_points_as_tori::generate::sampling::{center_and_scale, sample_mesh_uniform};
use ddg_points_as_tori::mc::{Band, McMesh, McOptions, MeshCheck, Resolution, marching_cubes};
use ddg_points_as_tori::spatial::{KdTree, KdTreeOptions};
use ddg_points_as_tori::truth::{MeshTruth, SignMode, TruthOptions, check_mesh};
use ddg_points_as_tori::{
    BuildOptions, EvalOptions, EvalOut, LambdaRule, Model, Neighborhood, NormalPolicy, PointCloud, PointsAsTori, SignRule, TorusField,
};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Instant;
use util::*;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let r = match args.first().map(String::as_str) {
        Some("prepare") => prepare(&args[1..]),
        Some("case") => case(&args[1..]),
        Some("continuity") => continuity(&args[1..]),
        Some("orient") => orient(&args[1..]),
        Some("selfint") => selfint(&args[1..]),
        Some("gpu") => gpu(&args[1..]),
        _ => Err("usage: eval prepare|case|continuity|orient ... (see the file header)".into()),
    };
    if let Err(e) = r {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

// ------------------------------------------------------------------ arguments

struct Args {
    pos: Vec<String>,
    opts: BTreeMap<String, Vec<String>>,
}

fn parse(a: &[String]) -> Args {
    let mut pos = Vec::new();
    let mut opts: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut i = 0;
    while i < a.len() {
        if let Some(k) = a[i].strip_prefix("--") {
            if i + 1 < a.len() && !a[i + 1].starts_with("--") {
                opts.entry(k.into()).or_default().push(a[i + 1].clone());
                i += 2;
            } else {
                opts.entry(k.into()).or_default().push(String::new());
                i += 1;
            }
        } else {
            pos.push(a[i].clone());
            i += 1;
        }
    }
    Args { pos, opts }
}

impl Args {
    fn get<T: std::str::FromStr>(&self, k: &str, d: T) -> R<T> {
        match self.opts.get(k).and_then(|v| v.last()) {
            None => Ok(d),
            Some(s) => s.parse().map_err(|_| format!("--{k} {s}: not a valid value")),
        }
    }
    fn flag(&self, k: &str) -> bool {
        self.opts.contains_key(k)
    }
    fn pos(&self, i: usize, what: &str) -> R<&str> {
        self.pos.get(i).map(String::as_str).ok_or_else(|| format!("missing argument <{what}>"))
    }
}

fn threads(a: &Args) -> R<Option<usize>> {
    let t: usize = a.get("threads", 0)?;
    Ok(if t == 0 { None } else { Some(t) })
}

fn name_seed(s: &str) -> u64 {
    // FNV-1a
    let mut h: u64 = 0xcbf29ce484222325;
    for b in s.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn crate_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn repo() -> PathBuf {
    crate_dir().join("../..")
}

// ------------------------------------------------------------------ prepare

struct Meta {
    kv: Vec<(String, String)>,
}

impl Meta {
    fn new() -> Meta {
        Meta { kv: Vec::new() }
    }
    fn set(&mut self, k: &str, v: impl ToString) -> &mut Self {
        self.kv.push((k.into(), v.to_string()));
        self
    }
    fn write(&self, dir: &Path) -> R<()> {
        let s: String = self.kv.iter().map(|(k, v)| format!("{k}={v}\n")).collect();
        std::fs::write(dir.join("meta.txt"), s).map_err(err)
    }
    fn read(dir: &Path) -> R<BTreeMap<String, String>> {
        let s = std::fs::read_to_string(dir.join("meta.txt")).map_err(|e| format!("{}/meta.txt: {e}", dir.display()))?;
        Ok(s.lines().filter_map(|l| l.split_once('=')).map(|(k, v)| (k.to_string(), v.to_string())).collect())
    }
}

fn write_case(dir: &Path, p: &[[f64; 3]], n: &[[f64; 3]], mesh: Option<(&[[f64; 3]], &[[u32; 3]])>, meta: &Meta) -> R<()> {
    std::fs::create_dir_all(dir).map_err(err)?;
    write3(&dir.join("points.npy"), p)?;
    write3(&dir.join("normals.npy"), n)?;
    if let Some((v, f)) = mesh {
        write3(&dir.join("mesh_vertices.npy"), v)?;
        write_faces(&dir.join("mesh_faces.npy"), f)?;
    }
    meta.write(dir)
}

/// Normalize (vertex mean, max |coordinate| = 1), validate as a closed outward manifold, sample.
fn mesh_case(dir: &Path, v: &[[f64; 3]], f: &[[u32; 3]], npts: usize, seed: u64, meta: &mut Meta) -> R<()> {
    let v = center_and_scale(v);
    let rep = check_mesh(&v, f).map_err(err)?;
    if !rep.is_closed() {
        return Err(format!("not closed/manifold/oriented: {rep:?}"));
    }
    if rep.degenerate_faces > 0 {
        return Err(format!("{} degenerate faces", rep.degenerate_faces));
    }
    if !(rep.signed_volume > 0.0) {
        return Err(format!("signed volume {} (inward or flat)", rep.signed_volume));
    }
    let tm = ddg_core::topology::TriMesh::new(v.len(), f.to_vec()).map_err(|e| format!("TriMesh: {e}"))?;
    let mut rng = Rng::new(seed);
    let (p, n) = sample_mesh_uniform(&v, f, npts, &mut rng).map_err(err)?;
    meta.set("truth", "mesh")
        .set("faces", f.len())
        .set("vertices", v.len())
        .set("volume", rep.signed_volume)
        .set("components", tm.n_components())
        .set("chi", tm.euler_characteristic())
        .set("seed", seed);
    write_case(dir, &p, &n, Some((&v, f)), meta)
}

fn prepare(a: &[String]) -> R<()> {
    let a = parse(a);
    let ds = a.pos(0, "dataset")?.to_string();
    let root = PathBuf::from(a.pos(1, "root")?);
    let npts: usize = a.get("points", 512)?;
    let th = threads(&a)?;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(th.unwrap_or(0)).build().map_err(err)?;
    let tag = if npts == 512 { ds.clone() } else { format!("{ds}{npts}") };
    let out = root.join(&tag);
    std::fs::create_dir_all(&out).map_err(err)?;
    let src = root.join("src");
    let mut log: Vec<String> = Vec::new();
    match ds.as_str() {
        "abc" => {
            // 20 validation shapes (minus exact training duplicates, which are kept but flagged),
            // then seen shapes from the training files in manifest order until 100 are valid.
            let man = std::fs::read_to_string(src.join("abc/manifest.tsv")).map_err(|e| format!("abc manifest: {e} (run scripts/pat_eval_data.py)"))?;
            let overlap = std::fs::read_to_string(src.join("overlap.json")).unwrap_or_default();
            let rows: Vec<Vec<String>> = man.lines().skip(1).map(|l| l.split('\t').map(String::from).collect()).collect();
            let res: Vec<(usize, String, R<()>)> = pool.install(|| {
                use rayon::prelude::*;
                rows.par_iter()
                    .enumerate()
                    .map(|(i, r)| {
                        let name = r[0].clone();
                        let d = src.join("abc").join(&name);
                        let go = || -> R<()> {
                            let v = read3(&d.join("vertices.npy"))?;
                            let f = read_faces(&d.join("faces.npy"))?;
                            let dup = overlap.contains(&format!("\"{name}\": {{\"exact\": [\""));
                            let training = if r[3] == "unseen" && dup { "duplicate-of-training" } else { r[3].as_str() };
                            let mut m = Meta::new();
                            m.set("dataset", "abc").set("name", &name).set("source", format!("{}:{}", r[1], r[2])).set("training", training);
                            mesh_case(&out.join(&name), &v, &f, npts, name_seed(&name), &mut m)
                        };
                        let r = go();
                        (i, name, r)
                    })
                    .collect()
            });
            let mut kept_seen = 0;
            let mut kept = 0;
            for (_, name, r) in &res {
                let seen = !name.contains("val");
                match r {
                    Ok(()) if seen && kept_seen >= 80 => {
                        let _ = std::fs::remove_dir_all(out.join(name));
                        log.push(format!("{name}\tsurplus"));
                    }
                    Ok(()) => {
                        kept += 1;
                        if seen {
                            kept_seen += 1;
                        }
                        log.push(format!("{name}\tkept"));
                    }
                    Err(e) => {
                        let _ = std::fs::remove_dir_all(out.join(name));
                        log.push(format!("{name}\trejected: {e}"));
                    }
                }
            }
            log.push(format!("# kept {kept} ({kept_seen} seen in training)"));
        }
        "blobs" => {
            let res: Vec<(u64, R<()>)> = pool.install(|| {
                use rayon::prelude::*;
                (0..100u64)
                    .into_par_iter()
                    .map(|i| {
                        let go = || -> R<()> {
                            let (m, p) = generate_blob(14_014, i, &BlobOptions::default()).map_err(err)?;
                            let name = format!("blob{i:03}");
                            let mut meta = Meta::new();
                            meta.set("dataset", "blobs").set("name", &name).set("training", "unseen").set("params", format!("{p:?}"));
                            mesh_case(&out.join(&name), &m.vertices, &m.faces, npts, name_seed(&name), &mut meta)
                        };
                        (i, go())
                    })
                    .collect()
            });
            for (i, r) in res {
                log.push(format!("blob{i:03}\t{}", r.err().unwrap_or_else(|| "kept".into())));
            }
        }
        "thingi" => prepare_thingi(&out, npts, &pool, &mut log)?,
        "analytic" | "analyticsad" => {
            // App. E: 70 instances; `analyticsad`: the 13 with saddle regions (NOTES-WP8)
            let sad = [9u64, 10, 11, 13, 22, 23, 45, 46, 47, 48, 49, 55, 59];
            let seeds: Vec<u64> = if ds == "analytic" { (0..70).collect() } else { sad.to_vec() };
            for s in seeds {
                let inst = Instance::from_seed(s).map_err(err)?;
                let name = inst.name();
                let (p, n) = inst.solid.sample_surface(npts, name_seed(&name)).map_err(err)?;
                let mut m = Meta::new();
                m.set("dataset", &ds).set("name", &name).set("truth", format!("analytic:{s}")).set("training", "unseen");
                m.set("kind", PAPER_KINDS[(s as usize / 5) % PAPER_KINDS.len()].name());
                write_case(&out.join(&name), &p, &n, None, &m)?;
                log.push(format!("{name}\tkept"));
            }
        }
        "saddle" => {
            for s in saddle_set().map_err(err)? {
                let (p, n) = s.solid.sample_surface(npts, name_seed(s.name)).map_err(err)?;
                let mut m = Meta::new();
                m.set("dataset", "saddle").set("name", s.name).set("truth", format!("saddle:{}", s.name)).set("training", "unseen");
                write_case(&out.join(s.name), &p, &n, None, &m)?;
                log.push(format!("{}\tkept", s.name));
            }
        }
        "ref" => {
            // the reference's data/ clouds with a truth mesh (WP5's extraction), resampled
            for name in ["block", "knot", "femur"] {
                let d = repo().join("runs_tmp/pat-wp5/meshes").join(name);
                let v = read3(&d.join("vertices.npy"))?;
                let f = read_faces(&d.join("faces.npy"))?;
                let mut m = Meta::new();
                m.set("dataset", "ref").set("name", name).set("training", "unseen");
                let r = mesh_case(&out.join(name), &v, &f, npts, name_seed(name), &mut m);
                log.push(format!("{name}\t{}", r.err().unwrap_or_else(|| "kept".into())));
            }
        }
        "reffull" => {
            // the 12 golden inputs as given (golden normalization), with the authors' golden grids
            let gold = repo().join("runs_tmp/pat-ref");
            for name in
                ["synth_sphere", "synth_torus", "ABC0_1", "ABC0_6", "ABC17_1", "ABC17_3", "block", "knot", "tooth", "femur", "snail", "holey-bunny"]
            {
                let inp = gold.join("inputs").join(name);
                let p = read3(&inp.join("points.npy"))?;
                let n = read3(&inp.join("normals.npy"))?;
                let mut m = Meta::new();
                m.set("dataset", "reffull").set("name", name).set("training", if name.starts_with("ABC17") { "validation-file" } else { "unseen" });
                let md = repo().join("runs_tmp/pat-wp5/meshes").join(name);
                let mesh = if md.join("faces.npy").is_file() { Some((read3(&md.join("vertices.npy"))?, read_faces(&md.join("faces.npy"))?)) } else { None };
                match name {
                    "synth_sphere" => {
                        m.set("truth", "sphere:0.5");
                    }
                    "synth_torus" => {
                        m.set("truth", "torus:0.35,0.15");
                    }
                    _ if mesh.is_some() => {
                        m.set("truth", "mesh");
                    }
                    _ => {
                        m.set("truth", "none");
                    }
                }
                let d = out.join(name);
                write_case(&d, &p, &n, mesh.as_ref().map(|(v, f)| (v.as_slice(), f.as_slice())), &m)?;
                // the golden outputs of the authors' program on exactly this cloud
                let g = gold.join("cpu_run1").join(name);
                let rd = d.join("ref");
                std::fs::create_dir_all(&rd).map_err(err)?;
                for f in ["grid_sdf_accel.npy", "grid_sdf_noaccel.npy", "net_coeffs.npy", "tori_centers.npy", "tori_axes_angles.npy", "tori_major_radius.npy", "tori_minor_radius.npy", "knn_indices.npy"] {
                    std::fs::copy(g.join(f), rd.join(f)).map_err(|e| format!("{}: {e}", g.join(f).display()))?;
                }
                log.push(format!("{name}\tkept"));
            }
        }
        "orient" => {
            // WP9a's 11 clouds (5 analytic solids, 6 golden clouds) at their sizes
            for s in [3u64, 12, 26, 40, 52] {
                let inst = Instance::from_seed(s).map_err(err)?;
                let name = inst.name();
                let (p, n) = inst.solid.sample_surface(2048, name_seed(&name)).map_err(err)?;
                let mut m = Meta::new();
                m.set("dataset", "orient").set("name", &name).set("truth", format!("analytic:{s}"));
                write_case(&out.join(&name), &p, &n, None, &m)?;
                log.push(format!("{name}\tkept"));
            }
            let gold = repo().join("runs_tmp/pat-ref/inputs");
            for name in ["synth_torus", "tooth", "knot", "block", "femur", "snail"] {
                let p = read3(&gold.join(name).join("points.npy"))?;
                let n = read3(&gold.join(name).join("normals.npy"))?;
                let mut m = Meta::new();
                m.set("dataset", "orient").set("name", name).set("truth", "none");
                write_case(&out.join(name), &p, &n, None, &m)?;
                log.push(format!("{name}\tkept"));
            }
        }
        _ => return Err(format!("unknown dataset {ds}")),
    }
    let txt: String = log.iter().map(|l| format!("{l}\n")).collect();
    std::fs::write(out.join("prepare.tsv"), &txt).map_err(err)?;
    print!("{txt}");
    Ok(())
}

/// Thingi10K selection (SPEC §16.2): closed, manifold, oriented, non-self-intersecting.
///
/// Candidates: our manifold conversion `data/thingi10k/converted/<id>.obj` (status `kept`, no
/// boundary) whose raw file the dataset's `geometry_data.csv` lists with 0 self-intersections,
/// vertex- and edge-manifold, oriented, no boundary, no degenerate or duplicated faces, between
/// 1,000 and 40,000 faces, and a licence that allows modified renders (CC BY, CC BY-SA, public
/// domain). Seeded shuffle; each candidate is then checked by `truth::check_mesh` (closed, outward,
/// no degenerate faces), `ddg_core::TriMesh::new` and our self-intersection test
/// ([`util::self_intersections`]: no edge crossing another face, no folded edge) until 45 pass.
fn prepare_thingi(out: &Path, npts: usize, pool: &rayon::ThreadPool, log: &mut Vec<String>) -> R<()> {
    let base = repo().join("data/thingi10k");
    let conv = std::fs::read_to_string(base.join("converted.tsv")).map_err(err)?;
    let mut kept_conv = std::collections::HashMap::new();
    for l in conv.lines().skip(1) {
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() > 8 && c[2] == "kept" && c[8] == "0" {
            kept_conv.insert(c[0].to_string(), c[7].parse::<usize>().unwrap_or(0));
        }
    }
    let geo = std::fs::read_to_string(base.join("metadata/geometry_data.csv")).map_err(err)?;
    let mut lic = std::collections::HashMap::new();
    let att = std::fs::read_to_string(base.join("ATTRIBUTION.tsv")).map_err(err)?;
    for l in att.lines().skip(1) {
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() >= 5 {
            lic.insert(c[0].to_string(), (c[4].to_string(), c[2].to_string(), c[3].to_string()));
        }
    }
    let ok_lic = |l: &str| {
        matches!(
            l,
            "Creative Commons - Attribution" | "Creative Commons - Attribution - Share Alike" | "Creative Commons - Public Domain Dedication" | "Public Domain"
        )
    };
    let mut cands: Vec<String> = Vec::new();
    for l in geo.lines().skip(1) {
        let c: Vec<&str> = l.split(',').collect();
        if c.len() < 14 {
            continue;
        }
        let num = |i: usize| c[i].parse::<f64>().unwrap_or(-1.0);
        let id = c[0].to_string();
        let Some(&faces) = kept_conv.get(&id) else { continue };
        let good = num(9) == 0.0 && num(11) == 1.0 && num(12) == 1.0 && num(13) == 1.0 && num(6) == 0.0 && num(3) == 0.0 && num(7) == 0.0;
        let l = lic.get(&id).map(|x| x.0.as_str()).unwrap_or("");
        if good && (1000..=40000).contains(&faces) && ok_lic(l) {
            cands.push(id);
        }
    }
    cands.sort();
    let mut rng = Rng::new(45);
    rng.shuffle(&mut cands);
    log.push(format!("# {} candidates after metadata filters", cands.len()));
    let mut kept = 0;
    for chunk in cands.chunks(24) {
        if kept >= 45 {
            break;
        }
        let res: Vec<(String, R<()>)> = pool.install(|| {
            use rayon::prelude::*;
            chunk
                .par_iter()
                .map(|id| {
                    let go = || -> R<()> {
                        let obj = ddg_core::io::read_obj(base.join("converted").join(format!("{id}.obj"))).map_err(err)?;
                        let f = obj.triangles().ok_or("not a triangle mesh")?;
                        let v: Vec<[f64; 3]> = obj.positions.clone();
                        let vn = center_and_scale(&v);
                        let si = self_intersections(&vn, &f);
                        if si.crossings > 0 || si.folds > 0 {
                            return Err(format!("self-intersecting: {} crossings, {} folds ({} pairs tested)", si.crossings, si.folds, si.pairs_tested));
                        }
                        let (l, nm, au) = lic.get(id).cloned().unwrap_or_default();
                        let mut m = Meta::new();
                        m.set("dataset", "thingi").set("name", format!("thingi{id}")).set("training", "unseen").set("licence", l).set("model", nm).set("author", au);
                        m.set("coplanar_close", si.coplanar_close).set("pairs_tested", si.pairs_tested);
                        mesh_case(&out.join(format!("thingi{id}")), &v, &f, npts, name_seed(id), &mut m)
                    };
                    (id.clone(), go())
                })
                .collect()
        });
        for (id, r) in res {
            match r {
                Ok(()) if kept < 45 => {
                    kept += 1;
                    log.push(format!("thingi{id}\tkept"));
                }
                Ok(()) => {
                    let _ = std::fs::remove_dir_all(out.join(format!("thingi{id}")));
                    log.push(format!("thingi{id}\tsurplus"));
                }
                Err(e) => {
                    let _ = std::fs::remove_dir_all(out.join(format!("thingi{id}")));
                    log.push(format!("thingi{id}\trejected: {e}"));
                }
            }
        }
    }
    log.push(format!("# kept {kept}"));
    Ok(())
}

// ------------------------------------------------------------------ truth

enum Truth {
    Mesh(MeshTruth, usize, i64),
    Solid(Solid),
    Sphere(f64),
    Torus(f64, f64),
    None,
}

fn load_truth(dir: &Path, meta: &BTreeMap<String, String>, th: Option<usize>, sign: SignMode) -> R<Truth> {
    let t = meta.get("truth").map(String::as_str).unwrap_or("none");
    Ok(if t == "mesh" {
        let v = read3(&dir.join("mesh_vertices.npy"))?;
        let f = read_faces(&dir.join("mesh_faces.npy"))?;
        let (mt, _) = MeshTruth::new(&v, &f, &TruthOptions { threads: th, sign, ..Default::default() }).map_err(err)?;
        let tm = ddg_core::topology::TriMesh::new(v.len(), f).map_err(|e| format!("TriMesh: {e}"))?;
        Truth::Mesh(mt, tm.n_components(), tm.euler_characteristic())
    } else if let Some(s) = t.strip_prefix("analytic:") {
        Truth::Solid(Instance::from_seed(s.parse().map_err(err)?).map_err(err)?.solid)
    } else if let Some(s) = t.strip_prefix("saddle:") {
        Truth::Solid(saddle_set().map_err(err)?.into_iter().find(|x| x.name == s).ok_or("unknown saddle shape")?.solid)
    } else if let Some(s) = t.strip_prefix("sphere:") {
        Truth::Sphere(s.parse().map_err(err)?)
    } else if let Some(s) = t.strip_prefix("torus:") {
        let (a, b) = s.split_once(',').ok_or("torus:R,r")?;
        Truth::Torus(a.parse().map_err(err)?, b.parse().map_err(err)?)
    } else {
        Truth::None
    })
}

fn eval_truth(t: &Truth, xs: &[[f64; 3]], th: Option<usize>) -> R<Option<Vec<f64>>> {
    use rayon::prelude::*;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(th.unwrap_or(0)).build().map_err(err)?;
    Ok(match t {
        Truth::Mesh(m, _, _) => {
            let mut s = vec![0.0; xs.len()];
            m.evaluate(xs, &mut s, None).map_err(err)?;
            Some(s)
        }
        Truth::Solid(s) => Some(pool.install(|| xs.par_iter().map(|&x| s.eval(x).0).collect())),
        Truth::Sphere(r) => Some(xs.iter().map(|x| (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt() - r).collect()),
        Truth::Torus(big, small) => Some(xs.iter().map(|x| ((x[0] * x[0] + x[1] * x[1]).sqrt() - big).hypot(x[2]) - small).collect()),
        Truth::None => None,
    })
}

fn load_cloud(dir: &Path) -> R<PointCloud> {
    let p = read3(&dir.join("points.npy"))?;
    let n = read3(&dir.join("normals.npy"))?;
    Ok(PointCloud::new(p, n, NormalPolicy::RequireUnit { tol: 1e-6 }).map_err(err)?.0)
}

fn pat_eval(field: &TorusField, xs: &[[f64; 3]]) -> R<(Vec<f64>, f64)> {
    let mut phi = vec![0.0; xs.len()];
    let rep = field.evaluate(xs, EvalOut { distance: Some(&mut phi), ..Default::default() }).map_err(err)?;
    Ok((phi, rep.seconds))
}

fn slices(n: usize, phi: &[f64]) -> Vec<f64> {
    // the three mid-planes (x, y, z index n/2)
    let c = n / 2;
    let mut s = Vec::with_capacity(3 * n * n);
    for a in 0..n {
        for b in 0..n {
            s.push(phi[gi(n, c, a, b)]);
        }
    }
    for a in 0..n {
        for b in 0..n {
            s.push(phi[gi(n, a, c, b)]);
        }
    }
    for a in 0..n {
        for b in 0..n {
            s.push(phi[gi(n, a, b, c)]);
        }
    }
    s
}

// ------------------------------------------------------------------ case

fn case(a: &[String]) -> R<()> {
    let a = parse(a);
    let dir = PathBuf::from(a.pos(0, "case-dir")?);
    let th = threads(&a)?;
    let meta = Meta::read(&dir)?;
    let cloud = load_cloud(&dir)?;
    let n = 64;
    let xs = grid(n);
    let t0 = Instant::now();
    let truth = load_truth(&dir, &meta, th, SignMode::Auto)?;
    let tv = eval_truth(&truth, &xs, th)?;
    let truth_s = t0.elapsed().as_secs_f64();
    let model = Model::reference().map_err(err)?;
    let mut methods: Vec<(String, Vec<f64>, f64)> = Vec::new();
    let mut extra = Json::new();
    extra.num("truth_seconds", truth_s);
    // PAT: one build per sign rule, then several evaluation rules on the same tori
    let mut tori_ref = None;
    let mut build_seconds = 0.0;
    for (rule, tag) in [(SignRule::Reference, ""), (SignRule::Geometric, "-geo")] {
        let mut opts = BuildOptions::default().with_threads(th);
        opts.fit.sign_rule = rule;
        let (field, rep) = PointsAsTori::build(&cloud, &model, &opts).map_err(err)?;
        extra.int(&format!("saddles{tag}"), rep.fit.saddle as i64).int(&format!("zero_curvature{tag}"), rep.fit.zero_curvature as i64);
        if tag.is_empty() {
            build_seconds = rep.seconds;
            extra.num("pat_build_seconds", rep.seconds).num("pat_predict_seconds", rep.predict.forward_seconds + rep.predict.features_seconds);
        }
        let (phi, s) = pat_eval(&field, &xs)?;
        methods.push((format!("pat{tag}"), phi, s));
        let all = TorusField::new(&cloud, field.tori().clone(), EvalOptions { neighborhood: Neighborhood::All, threads: th, ..Default::default() }).map_err(err)?;
        let (phi, s) = pat_eval(&all, &xs)?;
        methods.push((format!("pat{tag}-all"), phi, s));
        if tag.is_empty() {
            let (po, _) = TorusField::paper_options(&cloud, th).map_err(err)?;
            let paper = TorusField::new(&cloud, field.tori().clone(), po).map_err(err)?;
            let (phi, s) = pat_eval(&paper, &xs)?;
            methods.push(("pat-paper".into(), phi, s));
            let pool = rayon::ThreadPoolBuilder::new().num_threads(th.unwrap_or(0)).build().map_err(err)?;
            let tf = TaperField { field: &field, k: 32, c: 64.0 };
            let mut phi = vec![0.0; xs.len()];
            let t = Instant::now();
            pool.install(|| tf.evaluate(&xs, &mut phi)).map_err(err)?;
            methods.push(("pat-taper".into(), phi, t.elapsed().as_secs_f64()));
            if a.flag("sweep") {
                for (lab, eo) in [
                    ("pat-c32", EvalOptions { lambda: LambdaRule::PerQuery { c: 32.0 }, threads: th, ..Default::default() }),
                    ("pat-c128", EvalOptions { lambda: LambdaRule::PerQuery { c: 128.0 }, threads: th, ..Default::default() }),
                    ("pat-k16", EvalOptions { neighborhood: Neighborhood::Knn { k: 16 }, threads: th, ..Default::default() }),
                    ("pat-k64", EvalOptions { neighborhood: Neighborhood::Knn { k: 64 }, threads: th, ..Default::default() }),
                    ("pat-all-c32", EvalOptions { neighborhood: Neighborhood::All, lambda: LambdaRule::PerQuery { c: 32.0 }, threads: th, ..Default::default() }),
                    ("pat-all-c128", EvalOptions { neighborhood: Neighborhood::All, lambda: LambdaRule::PerQuery { c: 128.0 }, threads: th, ..Default::default() }),
                ] {
                    let f = TorusField::new(&cloud, field.tori().clone(), eo).map_err(err)?;
                    let (phi, s) = pat_eval(&f, &xs)?;
                    methods.push((lab.into(), phi, s));
                }
            }
            tori_ref = Some(field);
        }
    }
    // retrained weights (WP10): --retrained label=path[:geometric]
    for spec in a.opts.get("retrained").cloned().unwrap_or_default() {
        let (label, rest) = spec.split_once('=').ok_or("--retrained label=path[:geometric]")?;
        let (path, geo) = match rest.strip_suffix(":geometric") {
            Some(p) => (p, true),
            None => (rest, false),
        };
        let bytes = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
        let m = Model::from_safetensors(&bytes).map_err(err)?;
        let mut opts = BuildOptions::default().with_threads(th);
        if geo {
            opts.fit.sign_rule = SignRule::Geometric;
        }
        let (field, _) = PointsAsTori::build(&cloud, &m, &opts).map_err(err)?;
        let (phi, s) = pat_eval(&field, &xs)?;
        methods.push((label.to_string(), phi, s));
    }
    // baselines (Barill areas, the PAT neighbor rule: 32-NN, λ = 64/σ_x)
    let t = Instant::now();
    let (areas, arep) = point_areas(&cloud, &PointAreaOptions { threads: th, ..Default::default() }).map_err(err)?;
    let area_s = t.elapsed().as_secs_f64();
    extra.int("areas_clipped", arep.clipped as i64).num("areas_total", arep.total_area);
    let kf = KernelField::new(&cloud, areas.clone(), KernelOptions { threads: th, shc_zero_value: Some(0.0), ..Default::default() }).map_err(err)?;
    for (lab, m) in [("shc", KernelMethod::Shc), ("sspd", KernelMethod::Sspd), ("lse", KernelMethod::Lse)] {
        let mut v = vec![0.0; xs.len()];
        let rep = kf.evaluate(m, &xs, &mut v, None).map_err(err)?;
        if lab == "shc" {
            extra.int("shc_zero_w", rep.degenerate as i64);
        }
        methods.push((lab.into(), v, rep.seconds));
    }
    let tree = KdTree::new(cloud.positions(), KdTreeOptions::default()).map_err(err)?;
    let t = Instant::now();
    let nearest: Vec<f64> = {
        use rayon::prelude::*;
        let pool = rayon::ThreadPoolBuilder::new().num_threads(th.unwrap_or(0)).build().map_err(err)?;
        pool.install(|| {
            xs.par_chunks(1024)
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
        })
    };
    let near_s = t.elapsed().as_secs_f64();
    let mut gwn = vec![0.0; xs.len()];
    let wr = winding_numbers(&cloud, &areas, &xs, &mut gwn, th).map_err(err)?;
    let mut rwn = vec![0.0; xs.len()];
    let rr = regularized_winding_numbers(&cloud, &areas, 0.02, &xs, &mut rwn, th).map_err(err)?;
    let lse = methods.iter().find(|m| m.0 == "lse").map(|m| (m.1.clone(), m.2)).ok_or("lse missing")?;
    let signed = |w: &[f64], d: &[f64]| -> Vec<f64> { w.iter().zip(d).map(|(&w, &d)| winding_signed_distance(d, w, 0.5)).collect() };
    methods.push(("gwn-nearest".into(), signed(&gwn, &nearest), wr.seconds + near_s + area_s));
    methods.push(("rwn-nearest".into(), signed(&rwn, &nearest), rr.seconds + near_s + area_s));
    methods.push(("lse-gwn".into(), signed(&gwn, &lse.0), wr.seconds + lse.1 + area_s));
    methods.retain(|m| m.0 != "lse");
    // SPSR mesh SDF (scripts/pat_spsr.py output in <case>/spsr)
    let sd = dir.join("spsr");
    if sd.join("faces.npy").is_file() {
        let v = read3(&sd.join("vertices.npy"))?;
        let f = read_faces(&sd.join("faces.npy"))?;
        if f.is_empty() {
            extra.str("spsr", "empty mesh");
        } else {
            let t = Instant::now();
            let (mt, rep) = MeshTruth::new(&v, &f, &TruthOptions { threads: th, require_closed: false, ..Default::default() }).map_err(err)?;
            let mut s = vec![0.0; xs.len()];
            mt.evaluate(&xs, &mut s, None).map_err(err)?;
            extra.int("spsr_faces", rep.faces as i64).int("spsr_boundary", rep.boundary_edges as i64).int("spsr_nonmanifold", rep.nonmanifold_edges as i64).num("spsr_volume", rep.signed_volume);
            methods.push(("spsr".into(), s, t.elapsed().as_secs_f64()));
        }
    }
    // the authors' program (scripts/pat_eval_ref.py output in <case>/ref)
    let rd = dir.join("ref");
    let mut cmp = Json::new();
    let mut have_ref = false;
    for (file, lab) in [("grid_sdf_accel.npy", "authors"), ("grid_sdf_noaccel.npy", "authors-all")] {
        if rd.join(file).is_file() {
            let (s, v) = read_flat(&rd.join(file))?;
            if s.iter().product::<usize>() != xs.len() {
                return Err(format!("{}: shape {s:?}", rd.join(file).display()));
            }
            methods.push((lab.into(), v, f64::NAN));
            have_ref = true;
        }
    }
    if have_ref {
        let field = tori_ref.as_ref().ok_or("no reference field")?;
        black_box(&rd, field, &methods, &mut cmp)?;
    }
    // metrics
    let mut mj = Json::new();
    let mut sl: Vec<f64> = Vec::new();
    let mut names = Vec::new();
    if let Some(tv) = &tv {
        sl.extend(slices(n, tv));
        names.push(jstr("truth"));
    }
    for (lab, phi, s) in &methods {
        let mut j = Json::new();
        if let Some(tv) = &tv {
            let m = metrics(n, phi, tv);
            j.raw("m", &metrics_json(&m));
        }
        let jp = jumps(n, phi, 0);
        j.num("seconds", *s).num("jump_max", jp.max).num("jump_p999", jp.p999).int("jump_over5", jp.over5 as i64);
        mj.raw(lab, &j.finish());
        sl.extend(slices(n, phi));
        names.push(jstr(lab));
    }
    write_npy_f32(&dir.join("slices.npy"), &[names.len(), 3, n, n], &sl)?;
    let mut out = Json::new();
    for (k, v) in &meta {
        out.str(k, v);
    }
    out.int("n", cloud.len() as i64)
        .str("threads", &format!("{th:?}"))
        .num("pat_build_s", build_seconds)
        .raw("extra", &extra.finish())
        .raw("methods", &mj.finish())
        .raw("slice_names", &jarr(&names));
    if have_ref {
        out.raw("black_box", &cmp.finish());
    }
    std::fs::write(dir.join("result.json"), out.finish()).map_err(err)?;
    println!("{} done in {:.1}s", dir.display(), t0.elapsed().as_secs_f64());
    Ok(())
}

/// Ours vs the authors' program on the same cloud: coefficients, tori, SDF.
fn black_box(rd: &Path, field: &TorusField, methods: &[(String, Vec<f64>, f64)], j: &mut Json) -> R<()> {
    let tori = field.tori();
    let n = tori.tori.len();
    let stat = |mut v: Vec<f64>| -> String {
        let mx = v.iter().copied().fold(0.0, f64::max);
        let med = quantile(&mut v, 0.5);
        let p999 = quantile(&mut v, 0.999);
        Json::new().num("max", mx).num("median", med).num("p999", p999).finish()
    };
    if rd.join("net_coeffs.npy").is_file() {
        let (s, c) = read_flat(&rd.join("net_coeffs.npy"))?;
        if s == [n, 6] {
            let d: Vec<f64> = (0..n * 6).map(|i| (tori.coefficients[i / 6][i % 6] as f64 - c[i]).abs() / (1.0 + c[i].abs())).collect();
            j.raw("coeff_rel", &stat(d));
        }
    }
    if rd.join("knn_indices.npy").is_file() {
        // count rows whose 64-NN set differs (as a set) from ours
        let (s, k) = read_flat(&rd.join("knn_indices.npy"))?;
        if s.len() == 2 && s[0] == n {
            let kk = s[1];
            let mut buf = Vec::new();
            let mut diff = 0;
            for i in 0..n {
                field.tree().knn(field.positions()[i], kk, &mut buf);
                let mut a: Vec<u32> = buf.iter().map(|b| b.index).collect();
                let mut b: Vec<u32> = k[i * kk..(i + 1) * kk].iter().map(|&x| x as u32).collect();
                a.sort_unstable();
                b.sort_unstable();
                diff += (a != b) as usize;
            }
            j.int("knn_rows_differ", diff as i64);
        }
    }
    if rd.join("tori_centers.npy").is_file() {
        let c = read3(&rd.join("tori_centers.npy"))?;
        let (_, ang) = read_flat(&rd.join("tori_axes_angles.npy"))?;
        let (_, big) = read_flat(&rd.join("tori_major_radius.npy"))?;
        let (_, small) = read_flat(&rd.join("tori_minor_radius.npy"))?;
        if c.len() == n {
            let mut dc = Vec::new();
            let mut dr = Vec::new();
            let mut da = Vec::new();
            let mut sign_flips = 0;
            let mut flat = 0;
            for i in 0..n {
                let t = &tori.tori[i];
                if small[i].signum() != t.minor.signum() {
                    sign_flips += 1;
                }
                // near-flat tori (radii ≥ 1e3) are compared only through the SDF (pat-ref README)
                if big[i].abs() > 1e3 || small[i].abs() > 1e3 {
                    flat += 1;
                    continue;
                }
                let scale = 1.0 + big[i].abs() + small[i].abs();
                dc.push((0..3).map(|k| (t.center[k] - c[i][k]).abs()).fold(0.0, f64::max) / scale);
                dr.push(((t.major - big[i]).abs().max((t.minor - small[i]).abs())) / scale);
                let u = Torus::axis_from_angles(ang[2 * i], ang[2 * i + 1]);
                let cr = [t.axis[1] * u[2] - t.axis[2] * u[1], t.axis[2] * u[0] - t.axis[0] * u[2], t.axis[0] * u[1] - t.axis[1] * u[0]];
                da.push((cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt());
            }
            j.raw("tori_center_rel", &stat(dc)).raw("tori_radius_rel", &stat(dr)).raw("tori_axis_sin", &stat(da)).int("tori_sign_flips", sign_flips).int("tori_flat_skipped", flat);
        }
    }
    let get = |l: &str| methods.iter().find(|m| m.0 == l).map(|m| &m.1);
    for (ours, theirs, key) in [("pat", "authors", "sdf_accel"), ("pat-all", "authors-all", "sdf_all")] {
        if let (Some(a), Some(b)) = (get(ours), get(theirs)) {
            j.raw(key, &stat(a.iter().zip(b).map(|(x, y)| (x - y).abs()).collect()));
        }
    }
    if let Ok(t) = std::fs::read_to_string(rd.join("timing.json")) {
        j.raw("authors_timing", t.trim());
    }
    Ok(())
}

use ddg_points_as_tori::Torus;

// ------------------------------------------------------------------ continuity

/// Connected components of an MC mesh (by shared vertices) with their areas and the smallest
/// |φ_true| over their vertices.
fn components(m: &McMesh, truth: Option<&[f64]>) -> Vec<(f64, f64, usize)> {
    let nv = m.positions.len();
    let mut parent: Vec<u32> = (0..nv as u32).collect();
    fn find(p: &mut [u32], mut x: u32) -> u32 {
        while p[x as usize] != x {
            p[x as usize] = p[p[x as usize] as usize];
            x = p[x as usize];
        }
        x
    }
    for f in &m.faces {
        let a = find(&mut parent, f[0]);
        for &v in &f[1..] {
            let b = find(&mut parent, v);
            if a != b {
                parent[b as usize] = a;
            }
        }
    }
    let areas = m.areas();
    let mut acc: BTreeMap<u32, (f64, f64, usize)> = BTreeMap::new();
    for (fi, f) in m.faces.iter().enumerate() {
        let r = find(&mut parent, f[0]);
        let e = acc.entry(r).or_insert((0.0, f64::INFINITY, 0));
        e.0 += areas[fi];
        e.2 += 1;
    }
    if let Some(t) = truth {
        for v in 0..nv {
            let r = find(&mut parent, v as u32);
            if let Some(e) = acc.get_mut(&r) {
                e.1 = e.1.min(t[v].abs());
            }
        }
    }
    let mut v: Vec<(f64, f64, usize)> = acc.into_values().collect();
    v.sort_by(|a, b| b.0.total_cmp(&a.0));
    v
}

enum Rule<'a> {
    Pat(&'a TorusField),
    Taper(TaperField<'a>),
}

impl Rule<'_> {
    fn field(&self) -> &dyn ImplicitField {
        match self {
            Rule::Pat(f) => *f,
            Rule::Taper(t) => t,
        }
    }
}

fn continuity(a: &[String]) -> R<()> {
    let a = parse(a);
    let dir = PathBuf::from(a.pos(0, "case-dir")?);
    let th = threads(&a)?;
    let res: Vec<usize> = a.get::<String>("res", "64,128,256".into())?.split(',').map(|s| s.parse().map_err(err)).collect::<R<_>>()?;
    let top: usize = a.get("top", 40)?;
    let meta = Meta::read(&dir)?;
    let cloud = load_cloud(&dir)?;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(th.unwrap_or(0)).build().map_err(err)?;
    let truth = load_truth(&dir, &meta, th, SignMode::Fast { beta: 2.0 })?;
    let model = Model::reference().map_err(err)?;
    let t0 = Instant::now();
    let (knn, _) = PointsAsTori::build(&cloud, &model, &BuildOptions::default().with_threads(th)).map_err(err)?;
    let (po, r_eval) = TorusField::paper_options(&cloud, th).map_err(err)?;
    let paper = TorusField::new(&cloud, knn.tori().clone(), po).map_err(err)?;
    let all = TorusField::new(&cloud, knn.tori().clone(), EvalOptions { neighborhood: Neighborhood::All, threads: th, ..Default::default() }).map_err(err)?;
    let rules: Vec<(&str, Rule)> = vec![
        ("knn32", Rule::Pat(&knn)),
        ("paper", Rule::Pat(&paper)),
        ("all", Rule::Pat(&all)),
        ("taper", Rule::Taper(TaperField { field: &knn, k: 32, c: 64.0 })),
    ];
    let mut out = Json::new();
    for (k, v) in &meta {
        out.str(k, v);
    }
    out.int("n", cloud.len() as i64).num("r_eval", r_eval);
    if let Truth::Mesh(_, c, chi) = &truth {
        out.int("truth_components", *c as i64).int("truth_chi", *chi);
    }
    let mut per_res = Json::new();
    for &n in &res {
        let xs = grid(n);
        let tv = eval_truth(&truth, &xs, th)?;
        let mut rj = Json::new();
        // truth topology at this resolution (MC of the truth field)
        let mc_opts = McOptions { bbox: Some([[-1.1; 3], [1.1; 3]]), resolution: Resolution::Cells([n; 3]), band: Band::Dense, check_mesh: true, threads: th, ..Default::default() };
        let truth_mc = match &truth {
            Truth::Mesh(m, _, _) => Some(marching_cubes(&TruthField(m), &mc_opts).map_err(err)?),
            Truth::Solid(s) => Some(marching_cubes(s, &mc_opts).map_err(err)?),
            _ => None,
        };
        if let Some((m, rep)) = &truth_mc {
            if let MeshCheck::Manifold { components, euler_characteristic, .. } = rep.mesh {
                rj.int("truth_mc_components", components as i64).int("truth_mc_chi", euler_characteristic);
            }
            rj.int("truth_mc_triangles", m.faces.len() as i64);
        }
        for (lab, rule) in &rules {
            if *lab == "all" && n > 128 && cloud.len() > 4096 {
                continue;
            }
            let t = Instant::now();
            let mut phi = vec![0.0; xs.len()];
            match rule {
                Rule::Pat(f) => {
                    f.evaluate(&xs, EvalOut { distance: Some(&mut phi), ..Default::default() }).map_err(err)?;
                }
                Rule::Taper(tf) => pool.install(|| tf.evaluate(&xs, &mut phi)).map_err(err)?,
            }
            let eval_s = t.elapsed().as_secs_f64();
            let jp = jumps(n, &phi, top);
            let mut j = Json::new();
            j.raw("jumps", &jumps_json(&jp)).num("eval_seconds", eval_s);
            if let Some(tv) = &tv {
                j.raw("m", &metrics_json(&metrics(n, &phi, tv)));
                // where the large jumps are: |φ_true| at the edges with slope > 2
                let h = 2.0 / (n - 1) as f64;
                let mut d: Vec<f64> = Vec::new();
                let all_e = jumps_over(n, &phi, 2.0);
                for &(ia, _) in &all_e {
                    d.push(tv[ia].abs());
                }
                j.int("over2_edges", all_e.len() as i64);
                if !d.is_empty() {
                    let far = d.iter().filter(|&&x| x > 0.2).count();
                    j.num("over2_truth_dist_median", quantile(&mut d, 0.5)).num("over2_far_frac", far as f64 / all_e.len() as f64).num("h", h);
                }
            }
            // zero level set
            let (mesh, rep) = marching_cubes(rule.field(), &mc_opts).map_err(err)?;
            let tvm = match &truth {
                Truth::None => None,
                _ => eval_truth(&truth, &mesh.positions, th)?,
            };
            let comps = components(&mesh, tvm.as_deref());
            let h = 2.2 / n as f64;
            let detached = comps.iter().filter(|c| c.1 > 2.0 * h).count();
            let total_area: f64 = comps.iter().map(|c| c.0).sum();
            let small = comps.iter().filter(|c| c.0 < 0.01 * total_area).count();
            j.int("mc_components", comps.len() as i64).int("mc_detached", detached as i64).int("mc_small", small as i64).int("mc_triangles", mesh.faces.len() as i64);
            if let MeshCheck::Manifold { euler_characteristic, boundary_edges, .. } = rep.mesh {
                j.int("mc_chi", euler_characteristic).int("mc_boundary", boundary_edges as i64);
            }
            // mechanism: the largest jumps of the kNN and radius rules
            if matches!(*lab, "knn32" | "paper") && n == *res.iter().max().unwrap_or(&n) {
                if let Rule::Pat(f) = rule {
                    j.raw("mechanism", &mechanism(f, &knn, n, &jp.top, &truth)?);
                }
            }
            rj.raw(lab, &j.finish());
            eprintln!("{} {lab} n={n}: max slope {:.2} ({:.1}s)", dir.display(), jp.max, t.elapsed().as_secs_f64());
        }
        per_res.raw(&n.to_string(), &rj.finish());
    }
    out.raw("res", &per_res.finish()).num("seconds", t0.elapsed().as_secs_f64());
    std::fs::write(dir.join("continuity.json"), out.finish()).map_err(err)?;
    Ok(())
}

fn jumps_over(n: usize, phi: &[f64], thr: f64) -> Vec<(usize, usize)> {
    let h = 2.0 / (n - 1) as f64;
    let mut v = Vec::new();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let a = gi(n, i, j, k);
                for (ax, b) in [(0, (i + 1 < n).then(|| gi(n, i + 1, j, k))), (1, (j + 1 < n).then(|| gi(n, i, j + 1, k))), (2, (k + 1 < n).then(|| gi(n, i, j, k + 1)))] {
                    if let Some(b) = b
                        && (phi[b] - phi[a]).abs() / h > thr
                    {
                        v.push((a, ax));
                    }
                }
            }
        }
    }
    v
}

/// For each of the largest grid-edge jumps: locate the discontinuity by bisection, and compare
/// the neighbor sets and weights on its two sides.
fn mechanism(f: &TorusField, knn: &TorusField, n: usize, top: &[(f64, usize, usize)], truth: &Truth) -> R<String> {
    let axs = axis(n);
    let at = |a: usize| {
        let (i, j, k) = (a / (n * n), (a / n) % n, a % n);
        [axs[i], axs[j], axs[k]]
    };
    let h = 2.0 / (n - 1) as f64;
    let mut rows = Vec::new();
    let (mut swaps, mut far, mut cont) = (0, 0, 0);
    let mut ratio = Vec::new();
    let mut tdist = Vec::new();
    let mut wfrac = Vec::new();
    for &(slope, a, ax) in top {
        let x0 = at(a);
        let phi = |x: [f64; 3]| f.distance(x).unwrap_or(f64::NAN);
        let (mut lo, mut hi) = (0.0f64, 1.0f64);
        let p = |s: f64| {
            let mut x = x0;
            x[ax] += s * h;
            x
        };
        let (mut flo, mut fhi) = (phi(p(lo)), phi(p(hi)));
        for _ in 0..50 {
            let m = 0.5 * (lo + hi);
            let fm = phi(p(m));
            if (fm - flo).abs() >= (fhi - fm).abs() {
                hi = m;
                fhi = fm;
            } else {
                lo = m;
                flo = fm;
            }
        }
        let jump = (fhi - flo).abs();
        let (xl, xh) = (p(lo), p(hi));
        if jump < 0.25 * slope * h {
            cont += 1; // the change is spread over the edge (steep but continuous)
        }
        let (bl, sl) = f.blend_at(xl, BlendWant::Distance).map_err(err)?;
        let (bh, sh) = f.blend_at(xh, BlendWant::Distance).map_err(err)?;
        let mut a1 = sl.clone();
        let mut b1 = sh.clone();
        a1.sort_unstable();
        b1.sort_unstable();
        let changed = a1 != b1;
        swaps += changed as usize;
        // far-field indicator ρ_min/σ (0: at a data point, 2: all neighbors equidistant)
        let rmin = sh.iter().map(|&i| dist(xh, f.positions()[i as usize])).fold(f64::INFINITY, f64::min);
        let rr = rmin / bh.sigma;
        ratio.push(rr);
        if rr > 1.0 {
            far += 1;
        }
        // weight share of the entering point(s)
        let lam = bh.lambda;
        let w = |i: u32, x: [f64; 3], s: f64| (-lam * (dist(x, f.positions()[i as usize]) - s)).exp();
        let wsum: f64 = sh.iter().map(|&i| w(i, xh, bh.shift)).sum();
        let entering: Vec<u32> = b1.iter().filter(|i| a1.binary_search(i).is_err()).copied().collect();
        let leaving: Vec<u32> = a1.iter().filter(|i| b1.binary_search(i).is_err()).copied().collect();
        let win: f64 = entering.iter().map(|&i| w(i, xh, bh.shift)).sum::<f64>() / wsum;
        wfrac.push(win);
        let td = eval_truth(truth, &[xh], Some(1))?.map(|v| v[0]).unwrap_or(f64::NAN);
        tdist.push(td.abs());
        let _ = knn;
        rows.push(
            Json::new()
                .num("slope", slope)
                .num("jump", jump)
                .bool("set_changed", changed)
                .int("set_lo", sl.len() as i64)
                .int("set_hi", sh.len() as i64)
                .int("entering", entering.len() as i64)
                .int("leaving", leaving.len() as i64)
                .num("w_entering_frac", win)
                .num("rho_min_over_sigma", rr)
                .num("truth", td)
                .num("phi_lo", bl.distance)
                .num("phi_hi", bh.distance)
                .finish(),
        );
    }
    let m = top.len().max(1) as f64;
    Ok(Json::new()
        .int("edges", top.len() as i64)
        .num("set_change_frac", swaps as f64 / m)
        .num("far_field_frac", far as f64 / m)
        .num("continuous_frac", cont as f64 / m)
        .num("rho_min_over_sigma_median", quantile(&mut ratio, 0.5))
        .num("truth_dist_median", quantile(&mut tdist, 0.5))
        .num("w_entering_median", quantile(&mut wfrac, 0.5))
        .raw("rows", &jarr(&rows[..rows.len().min(12)]))
        .finish())
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

// ------------------------------------------------------------------ orientation

fn orient(a: &[String]) -> R<()> {
    let a = parse(a);
    let dir = PathBuf::from(a.pos(0, "case-dir")?);
    let th = threads(&a)?;
    let frac: f64 = a.get("flip", 0.3)?;
    let seed: u64 = a.get("seed", 11)?;
    let cloud = load_cloud(&dir)?;
    let truth_n = cloud.normals().to_vec();
    let n = cloud.len();
    let mut idx: Vec<usize> = (0..n).collect();
    let mut rng = Rng::new(seed ^ name_seed(&dir.display().to_string()));
    rng.shuffle(&mut idx);
    let m = (frac * n as f64).round() as usize;
    let mut flipped = truth_n.clone();
    for &i in &idx[..m] {
        flipped[i] = [-flipped[i][0], -flipped[i][1], -flipped[i][2]];
    }
    let model = Model::reference().map_err(err)?;
    let mut out = Json::new();
    out.str("case", &dir.display().to_string()).int("n", n as i64).num("flip", frac).int("flipped", m as i64);
    for iters in [10usize, 30] {
        let opts = OrientOptions { max_iterations: iters, threads: th, ..Default::default() };
        let (res, rep) = orient_normals(cloud.positions(), &flipped, &model, &opts).map_err(err)?;
        let wrong = res.iter().zip(&truth_n).filter(|(a, b)| a[0] * b[0] + a[1] * b[1] + a[2] * b[2] < 0.0).count();
        let fl: Vec<String> = rep.flips_per_iteration.iter().map(|x| x.to_string()).collect();
        out.raw(
            &format!("it{iters}"),
            &Json::new()
                .int("wrong", wrong as i64)
                .num("residual", wrong as f64 / n as f64)
                .bool("converged", rep.converged)
                .int("negative_final", rep.negative_final as i64)
                .raw("flips", &jarr(&fl))
                .num("seconds", rep.seconds)
                .finish(),
        );
        eprintln!("{} it{iters}: residual {:.4}", dir.display(), wrong as f64 / n as f64);
    }
    std::fs::write(dir.join("orient.json"), out.finish()).map_err(err)?;
    Ok(())
}

// ------------------------------------------------------------------ self-intersection check

/// `eval selfint <mesh.obj>...`: the Thingi10K selection test on single meshes (to check it
/// against the dataset's own `num_self_intersections`).
fn selfint(a: &[String]) -> R<()> {
    for p in a {
        let obj = ddg_core::io::read_obj(p).map_err(err)?;
        let f = obj.triangles().ok_or("not a triangle mesh")?;
        let v = center_and_scale(&obj.positions);
        let s = self_intersections(&v, &f);
        println!("{p}\tfaces {}\tcrossings {}\tfolds {}\tcoplanar_close {}\tpairs {}", f.len(), s.crossings, s.folds, s.coplanar_close, s.pairs_tested);
    }
    Ok(())
}

// ------------------------------------------------------------------ GPU matmul precision

/// `eval gpu <case-dir>... [--threads T]` (build with `--features cuda`): the network on the GPU
/// with `GpuMatmul::Tf32` (default) and `Split3`, against the CPU network; the fit and the blend
/// run on the CPU in f64 in all three, so the differences are the network's. Writes
/// `<case>/gpu.json`: coefficient differences, φ differences on the 64³ grid (32-NN and all
/// points), and the grid MAE of each against the truth.
fn gpu(a: &[String]) -> R<()> {
    use ddg_points_as_tori::Device;
    use ddg_points_as_tori::network::GpuMatmul;
    let a = parse(a);
    let th = threads(&a)?;
    let model = Model::reference().map_err(err)?;
    let xs = grid(64);
    for d in &a.pos {
        let dir = PathBuf::from(d);
        let meta = Meta::read(&dir)?;
        let cloud = load_cloud(&dir)?;
        let tv = eval_truth(&load_truth(&dir, &meta, th, SignMode::Auto)?, &xs, th)?;
        let mut fields = Vec::new();
        for (lab, dev, mm) in [("cpu", Device::Cpu, GpuMatmul::Tf32), ("tf32", Device::Cuda { ordinal: 0 }, GpuMatmul::Tf32), ("split3", Device::Cuda { ordinal: 0 }, GpuMatmul::Split3)] {
            let mut o = BuildOptions::default().with_threads(th);
            o.predict.device = dev;
            o.predict.gpu_matmul = mm;
            let (f, rep) = PointsAsTori::build(&cloud, &model, &o).map_err(err)?;
            let all = TorusField::new(&cloud, f.tori().clone(), EvalOptions { neighborhood: Neighborhood::All, threads: th, ..Default::default() }).map_err(err)?;
            let (p, _) = pat_eval(&f, &xs)?;
            let (pa, _) = pat_eval(&all, &xs)?;
            fields.push((lab, f, p, pa, rep.seconds));
        }
        let mut j = Json::new();
        for (lab, f, p, pa, s) in &fields {
            let mut k = Json::new();
            k.num("build_s", *s);
            if let Some(tv) = &tv {
                k.num("mae", metrics(64, p, tv).mae).num("mae_all", metrics(64, pa, tv).mae);
            }
            if *lab != "cpu" {
                let (_, fc, pc, pac, _) = &fields[0];
                let c0 = &fc.tori().coefficients;
                let c1 = &f.tori().coefficients;
                let dc = c0.iter().zip(c1).flat_map(|(x, y)| (0..6).map(move |i| (x[i] as f64 - y[i] as f64).abs() / (1.0 + (x[i] as f64).abs()))).fold(0.0, f64::max);
                let mut dp: Vec<f64> = p.iter().zip(pc).map(|(x, y)| (x - y).abs()).collect();
                let mut dpa: Vec<f64> = pa.iter().zip(pac).map(|(x, y)| (x - y).abs()).collect();
                let mx = dp.iter().copied().fold(0.0, f64::max);
                let mxa = dpa.iter().copied().fold(0.0, f64::max);
                let sf = f.tori().tori.iter().zip(&fc.tori().tori).filter(|(x, y)| x.minor.signum() != y.minor.signum()).count();
                k.num("coeff_rel_max", dc).num("phi_max", mx).num("phi_median", quantile(&mut dp, 0.5)).num("phi_p999", quantile(&mut dp, 0.999)).num("phi_all_max", mxa).num("phi_all_p999", quantile(&mut dpa, 0.999)).int("minor_sign_flips", sf as i64);
            }
            j.raw(lab, &k.finish());
        }
        std::fs::write(dir.join("gpu.json"), j.finish()).map_err(err)?;
        println!("{d}: {}", j.finish());
    }
    Ok(())
}
