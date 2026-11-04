//! Determinism audit (SPEC.md §5.4, §18 WP15): every public CPU entry point returns the same bits
//! at 1, 3, 8 and 32 threads and from run to run.
//!
//! Covered here: build + evaluation (`determinism.rs` has the long version), level-set sampling,
//! morphology, offsets and Booleans through marching cubes, normal orientation, PCA normals, the
//! generators (blobs, scans, sampling, FPS, simplification, 2D), the analytic sampler, the
//! baselines (areas, SHC/SSPD/LSE, winding numbers, FWN, mesh truth), the 2D variant. Each
//! result is compared through its `Debug` text, which prints every f64 with its shortest
//! round-trip digits (so equal text is equal bits, `-0.0` and NaN included).

mod common;

use common::{ellipsoid_cloud, lcg};
use ddg_points_as_tori::apps::{Morph, MorphOptions, SampleOptions, morphology, sample_level_set};
use ddg_points_as_tori::apps::orient::{OrientOptions, orient_normals, pca_normals_random_sign};
use ddg_points_as_tori::baselines::{FastWinding, FwnOptions, KernelField, KernelMethod, KernelOptions, PointAreaOptions, point_areas, winding_numbers};
use ddg_points_as_tori::csg::Csg;
use ddg_points_as_tori::generate::blob::{BlobOptions, generate_blob};
use ddg_points_as_tori::generate::bvh::TriBvh;
use ddg_points_as_tori::generate::d2::{Scan2dOptions, farthest_point_sampling_2d, random_instance_2d};
use ddg_points_as_tori::generate::dataset::{CloudMode, Layout, generate_blob_records, sample_cloud};
use ddg_points_as_tori::generate::rng::Rng;
use ddg_points_as_tori::generate::sampling::{farthest_point_sampling, sample_mesh_uniform};
use ddg_points_as_tori::generate::scanner::{ScanOptions, simulate_scan};
use ddg_points_as_tori::generate::simplify::simplify_closed;
use ddg_points_as_tori::mc::{McOptions, Resolution, marching_cubes};
use ddg_points_as_tori::pool::make_pool;
use ddg_points_as_tori::truth::{MeshTruth, SignMode, TruthOptions};
use ddg_points_as_tori::analytic::Instance;
use ddg_points_as_tori::{EvalOut, Model, PointCloud, PointsAsTori, TorusField};
use std::f64::consts::PI;

const THREADS: [usize; 3] = [3, 8, 32];

/// `f(1)` twice, then `f(t)` for t in 3, 8, 32: all equal.
fn audit(name: &str, f: impl Fn(usize) -> String) {
    let base = f(1);
    assert!(!base.is_empty(), "{name}: empty result");
    assert!(base == f(1), "{name}: two runs on 1 thread differ");
    for t in THREADS {
        let r = f(t);
        if base != r {
            let at = base.bytes().zip(r.bytes()).position(|(a, b)| a != b).unwrap_or(base.len().min(r.len()));
            let lo = at.saturating_sub(40);
            panic!("{name}: {t} threads differ from 1 thread at byte {at}: ...{} | ...{}", &base[lo..(at + 40).min(base.len())], &r[lo..(at + 40).min(r.len())]);
        }
    }
    eprintln!("determinism: {name}: {} bytes of result identical at 1 (twice), 3, 8, 32 threads", base.len());
}

fn in_pool<R: Send>(t: usize, f: impl FnOnce() -> R + Send) -> R {
    make_pool(Some(t)).unwrap().install(f)
}

fn uv_sphere(nu: usize, nv: usize, r: f64) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut v = vec![[0.0, r, 0.0]];
    for j in 1..nv {
        let th = PI * j as f64 / nv as f64;
        for i in 0..nu {
            let ph = 2.0 * PI * i as f64 / nu as f64;
            v.push([r * th.sin() * ph.cos(), r * th.cos(), r * th.sin() * ph.sin()]);
        }
    }
    v.push([0.0, -r, 0.0]);
    let south = (v.len() - 1) as u32;
    let ring = |j: usize, i: usize| 1 + (j - 1) * nu + (i % nu);
    let mut f = Vec::new();
    for i in 0..nu {
        f.push([0, ring(1, i + 1) as u32, ring(1, i) as u32]);
        f.push([south, ring(nv - 1, i) as u32, ring(nv - 1, i + 1) as u32]);
    }
    for j in 1..nv - 1 {
        for i in 0..nu {
            let (a, b, c, d) = (ring(j, i) as u32, ring(j, i + 1) as u32, ring(j + 1, i + 1) as u32, ring(j + 1, i) as u32);
            f.push([a, b, c]);
            f.push([a, c, d]);
        }
    }
    (v, f)
}

fn build(cloud: &PointCloud, model: &Model, t: usize) -> TorusField {
    let mut b = common::ref_build().with_threads(Some(t));
    b.eval = TorusField::paper_options(cloud, Some(t)).unwrap().0;
    PointsAsTori::build(cloud, model, &b).unwrap().0
}

fn field_text(f: &TorusField) -> String {
    let qs: Vec<[f64; 3]> = {
        let mut s = 9;
        (0..400).map(|_| [1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8]).collect()
    };
    let (mut d, mut g, mut l) = (vec![0.0; qs.len()], vec![[0.0; 3]; qs.len()], vec![0.0; qs.len()]);
    f.evaluate(&qs, EvalOut { distance: Some(&mut d), gradient: Some(&mut g), laplacian: Some(&mut l) }).unwrap();
    format!("{:?}{:?}{:?}{:?}", f.tori(), d, g, l)
}

#[test]
fn build_and_evaluate() {
    let model = Model::reference().unwrap();
    let cloud = ellipsoid_cloud(800, 3, [0.6, 0.45, 0.3]);
    audit("build + evaluate (paper options)", |t| field_text(&build(&cloud, &model, t)));
}

#[test]
fn level_set_sampling_and_marching_cubes() {
    let model = Model::reference().unwrap();
    let cloud = ellipsoid_cloud(800, 3, [0.6, 0.45, 0.3]);
    let field = build(&cloud, &model, 8);
    let mc = |t| McOptions { resolution: Resolution::Longest(28), threads: Some(t), refine_steps: 2, normals: true, ..Default::default() };
    audit("sample_level_set", |t| {
        let o = SampleOptions { count: Some(1200), mc: mc(t), threads: Some(t), seed: 5, ..Default::default() };
        let (p, n, rep) = sample_level_set(&field, 0.02, &o).unwrap();
        format!("{p:?}{n:?}{:?}{:?}{:?}{:?}", rep.produced, rep.rejected, rep.residual_max, rep.residual_before_max)
    });
    audit("marching_cubes", |t| {
        let (m, r) = marching_cubes(&field, &mc(t)).unwrap();
        format!("{:?}{:?}", m, r.cells_visited)
    });
    // offset and Booleans of two fitted fields
    let shifted = {
        let p: Vec<[f64; 3]> = cloud.positions().iter().map(|p| [p[0] + 0.25, p[1], p[2]]).collect();
        PointCloud::new(p, cloud.normals().to_vec(), Default::default()).unwrap().0
    };
    audit("csg union/subtract/offset + marching_cubes", |t| {
        let a = Csg::field(build(&cloud, &model, t));
        let b = Csg::field(build(&shifted, &model, t));
        let u = a.union(b.offset(0.01)).subtract(Csg::field(build(&cloud, &model, t)).offset(-0.1));
        format!("{:?}", marching_cubes(&u, &mc(t)).unwrap().0)
    });
}

#[test]
fn morphology_and_orientation() {
    let model = Model::reference().unwrap();
    let cloud = ellipsoid_cloud(800, 3, [0.6, 0.45, 0.3]);
    let field = build(&cloud, &model, 8);
    audit("morphology (close)", |t| {
        let opts = MorphOptions {
            sample: SampleOptions { count: Some(600), mc: McOptions { resolution: Resolution::Longest(24), threads: Some(t), ..Default::default() }, threads: Some(t), seed: 2, ..Default::default() },
            build: common::ref_build().with_threads(Some(t)),
            paper_eval: true,
        };
        let (f, _, rep) = morphology(&field, &model, Morph::Close(0.03), &opts).unwrap();
        format!("{}{:?}", field_text(&f), rep.steps.iter().map(|s| (s.level, s.sample.produced, s.sample.rejected)).collect::<Vec<_>>())
    });
    // a cloud with 30 % of its normals flipped
    let mut s = 4;
    let normals: Vec<[f64; 3]> = cloud.normals().iter().map(|n| if lcg(&mut s) < 0.3 { [-n[0], -n[1], -n[2]] } else { *n }).collect();
    audit("orient_normals", |t| {
        let o = OrientOptions { max_iterations: 3, threads: Some(t), build: common::ref_build().with_threads(Some(t)), ..Default::default() };
        let (n, rep) = orient_normals(cloud.positions(), &normals, &model, &o).unwrap();
        format!("{:?}{:?}{:?}{:?}", n, rep.flips_per_iteration, rep.negative_per_iteration, rep.undecided)
    });
    audit("pca_normals_random_sign", |t| format!("{:?}", pca_normals_random_sign(cloud.positions(), 16, 5, Some(t)).unwrap()));
}

#[test]
fn baselines() {
    let cloud = ellipsoid_cloud(1000, 31, [0.8, 0.6, 0.5]);
    let (v, f) = uv_sphere(40, 20, 0.7);
    let mut s = 8;
    let xs: Vec<[f64; 3]> = (0..1500).map(|_| [2.0 * lcg(&mut s) - 1.0, 2.0 * lcg(&mut s) - 1.0, 2.0 * lcg(&mut s) - 1.0]).collect();
    audit("baselines (areas, SHC/SSPD/LSE, winding, FWN, mesh truth)", |t| {
        let th = Some(t);
        let (a, _) = point_areas(&cloud, &PointAreaOptions { threads: th, ..Default::default() }).unwrap();
        let mut out = format!("{a:?}");
        let k = KernelField::new(&cloud, a.clone(), KernelOptions { threads: th, ..Default::default() }).unwrap();
        for m in [KernelMethod::Shc, KernelMethod::Sspd, KernelMethod::Lse] {
            let (mut val, mut g) = (vec![0.0; xs.len()], vec![[0.0; 3]; xs.len()]);
            k.evaluate(m, &xs, &mut val, Some(&mut g)).unwrap();
            out += &format!("{val:?}{g:?}");
        }
        let mut w = vec![0.0; xs.len()];
        winding_numbers(&cloud, &a, &xs, &mut w, th).unwrap();
        out += &format!("{w:?}");
        FastWinding::from_cloud(&cloud, &a, &FwnOptions { threads: th, ..Default::default() }).unwrap().evaluate(&xs, &mut w).unwrap();
        out += &format!("{w:?}");
        for sign in [SignMode::Exact, SignMode::Fast { beta: 2.0 }] {
            let (tr, _) = MeshTruth::new(&v, &f, &TruthOptions { threads: th, sign, ..Default::default() }).unwrap();
            let (mut d, mut ww) = (vec![0.0; xs.len()], vec![0.0; xs.len()]);
            tr.evaluate(&xs, &mut d, Some(&mut ww)).unwrap();
            out += &format!("{d:?}{ww:?}");
        }
        out
    });
}

#[test]
fn generators() {
    let (v, f) = uv_sphere(30, 16, 0.8);
    audit("generate_blob + sample_cloud", |t| {
        in_pool(t, || {
            let (mesh, params) = generate_blob(7, 1, &BlobOptions::default()).unwrap();
            let (p, n) = sample_cloud(&mesh, &CloudMode::Scan(ScanOptions { spacing: (0.05, 0.1), ..Default::default() }), 800, 3).unwrap();
            format!("{:?}{:?}{:?}", params, p, n)
        })
    });
    audit("generate_blob_records", |t| {
        in_pool(t, || {
            let mode = CloudMode::Scan(ScanOptions { spacing: (0.05, 0.1), ..ScanOptions::default() });
            format!("{:?}", generate_blob_records(7, 0, 3, &mode, &Layout::default(), &BlobOptions::default()).unwrap())
        })
    });
    audit("simulate_scan, uniform sampling, FPS, simplification", |t| {
        in_pool(t, || {
            let bvh = TriBvh::new(&v, &f).unwrap();
            let scan = simulate_scan(&bvh, &ScanOptions { spacing: (0.05, 0.1), ..Default::default() }, 5).unwrap();
            let (p, n) = sample_mesh_uniform(&v, &f, 3000, &mut Rng::new(2)).unwrap();
            let fps = farthest_point_sampling(&p, 200, 0);
            let simp = simplify_closed(&v, &f, 120);
            format!("{scan:?}{p:?}{n:?}{fps:?}{simp:?}")
        })
    });
    audit("2D generators", |t| {
        in_pool(t, || {
            let inst = random_instance_2d(3, 4).unwrap();
            let cloud = inst.scan(&Scan2dOptions::default(), 6).unwrap();
            let fps = farthest_point_sampling_2d(&inst.sample_uniform(2000, &mut Rng::new(1)).unwrap().positions, 100, 0);
            format!("{cloud:?}{fps:?}")
        })
    });
    audit("analytic sampling", |t| {
        in_pool(t, || {
            let inst = Instance::from_seed(23).unwrap();
            format!("{:?}", inst.solid.sample_surface(2000, 99).unwrap())
        })
    });
}
