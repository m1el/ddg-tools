//! Generators (SPEC.md §10.3, §11): Perlin noise, blobs, simplification, ray casting, the
//! simulated scanner, farthest-point sampling, mesh sampling, the npz writer and the 2D
//! generators. Run with `--release`.

use ddg_points_as_tori::analytic::{Kind, Shape2d};
use ddg_points_as_tori::generate::blob::{BlobField, BlobOptions, BlobParams, generate_blob, generate_blob_mesh};
use ddg_points_as_tori::generate::bvh::TriBvh;
use ddg_points_as_tori::generate::d2::{Blob2dParams, Curve2d, Instance2d, Scan2dOptions, farthest_point_sampling_2d, noise_blob, random_instance_2d, random_polygon, simulate_scan_2d};
use ddg_points_as_tori::generate::dataset::{CloudMode, Layout, ShapeRecord, generate_blob_records, sample_cloud, write_training_npz};
use ddg_points_as_tori::generate::perlin::Perlin;
use ddg_points_as_tori::generate::rng::Rng;
use ddg_points_as_tori::generate::sampling::{center_and_scale, farthest_point_sampling, sample_mesh_uniform};
use ddg_points_as_tori::generate::scanner::{ScanOptions, fibonacci_sphere, simulate_scan};
use ddg_points_as_tori::generate::simplify::simplify_closed;
use ddg_points_as_tori::mc::{Band, McOptions, Resolution, marching_cubes};
use ddg_points_as_tori::npy::parse_npz_stored;
use ddg_points_as_tori::pool::make_pool;
use ddg_points_as_tori::truth::{MeshTruth, TruthOptions, check_mesh};
use std::f64::consts::PI;

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

fn ks_uniform(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len() as f64;
    v.iter().enumerate().map(|(i, &x)| ((i as f64 + 1.0) / n - x).abs().max((x - i as f64 / n).abs())).fold(0.0, f64::max)
}

#[test]
fn rng_is_deterministic_uniform_and_streams_differ() {
    let mut a = Rng::new(5);
    let mut b = Rng::new(5);
    assert!((0..100).all(|_| a.next_u64() == b.next_u64()));
    assert_ne!(Rng::for_item(1, 2).next_u64(), Rng::for_item(1, 3).next_u64());
    let mut r = Rng::new(1);
    let u: Vec<f64> = (0..50_000).map(|_| r.f64()).collect();
    assert!(ks_uniform(u) < 1.95 / (50_000f64).sqrt());
    let g: Vec<f64> = (0..50_000).map(|_| r.normal()).collect();
    let m = g.iter().sum::<f64>() / g.len() as f64;
    let var = g.iter().map(|x| (x - m).powi(2)).sum::<f64>() / g.len() as f64;
    assert!(m.abs() < 0.02 && (var - 1.0).abs() < 0.03);
    let rot = Rng::new(3).rotation();
    for i in 0..3 {
        for j in 0..3 {
            let d: f64 = (0..3).map(|k| rot[i][k] * rot[j][k]).sum();
            assert!((d - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
        }
    }
}

#[test]
fn perlin_noise_properties() {
    let p = Perlin::new(7);
    // zero at lattice points
    for (x, y, z) in [(0.0, 0.0, 0.0), (3.0, -2.0, 5.0), (10.0, 1.0, -7.0)] {
        assert!(p.noise3(x, y, z).abs() < 1e-12);
    }
    assert!(p.noise2(4.0, -3.0).abs() < 1e-12);
    let mut rng = Rng::new(1);
    let vals: Vec<f64> = (0..20_000).map(|_| p.noise3(rng.uniform(-20.0, 20.0), rng.uniform(-20.0, 20.0), rng.uniform(-20.0, 20.0))).collect();
    let m = vals.iter().sum::<f64>() / vals.len() as f64;
    let sd = (vals.iter().map(|x| (x - m).powi(2)).sum::<f64>() / vals.len() as f64).sqrt();
    assert!(vals.iter().all(|v| v.abs() <= 1.0), "noise left [-1, 1]");
    assert!(m.abs() < 0.02 && (0.15..0.45).contains(&sd), "mean {m} sd {sd}");
    // continuity: a small step changes the value a little
    let lip = (0..2000)
        .map(|_| {
            let x = [rng.uniform(-5.0, 5.0), rng.uniform(-5.0, 5.0), rng.uniform(-5.0, 5.0)];
            (p.noise3(x[0] + 1e-6, x[1], x[2]) - p.noise3(x[0], x[1], x[2])).abs() / 1e-6
        })
        .fold(0.0, f64::max);
    assert!(lip < 8.0, "gradient {lip}");
    // bases differ; same base is identical; the base is taken mod 256
    assert_ne!(Perlin::new(1).noise3(0.3, 0.4, 0.5), Perlin::new(2).noise3(0.3, 0.4, 0.5));
    assert_eq!(Perlin::new(1).noise3(0.3, 0.4, 0.5), Perlin::new(257).noise3(0.3, 0.4, 0.5));
    // fbm is normalized
    let f = (0..5000).map(|_| p.fbm3([rng.uniform(-9.0, 9.0), rng.uniform(-9.0, 9.0), rng.uniform(-9.0, 9.0)], 4, 0.5, 2.0).abs()).fold(0.0, f64::max);
    assert!(f <= 1.0);
}

#[test]
fn blobs_are_watertight_deterministic_and_within_budget() {
    let opts = BlobOptions::default();
    for index in 0..12u64 {
        let (m, params) = generate_blob(42, index, &opts).unwrap();
        assert!(m.faces.len() <= 2048 && m.faces.len() > 100, "{} faces", m.faces.len());
        let rep = check_mesh(&m.vertices, &m.faces).unwrap();
        assert!(rep.is_closed() && rep.degenerate_faces == 0 && rep.signed_volume > 0.0, "blob {index} {params:?}: {rep:?}");
        // centered and scaled into [-1,1]^3 touching the box
        let mx = m.vertices.iter().flat_map(|v| v.iter()).fold(0.0f64, |a, &c| a.max(c.abs()));
        assert!((mx - 1.0).abs() < 1e-12);
        let (m2, p2) = generate_blob(42, index, &opts).unwrap();
        assert_eq!((m.clone(), params), (m2, p2));
    }
    let a = generate_blob(42, 0, &opts).unwrap().0;
    let b = generate_blob(43, 0, &opts).unwrap().0;
    assert_ne!(a, b);
}

#[test]
fn blob_mc_matches_the_implicit_function() {
    let params = BlobParams { alpha: 0.4, scale: 1.3, octaves: 3, persistence: 0.5, lacunarity: 2.0, base_radius: 1.0, noise_base: 9, grid: 64 };
    let m = generate_blob_mesh(&params, &BlobOptions { max_faces: usize::MAX, ..BlobOptions::default() }).unwrap();
    let f = BlobField::new(params.clone());
    // vertices are on F = 0 up to the linear interpolation error of a 64^3 grid (median test:
    // the normalized distance is large only where ∇F is small)
    let h = 2.0 * 1.2 * 1.4 / 63.0;
    let mut errs: Vec<f64> = m
        .vertices
        .iter()
        .map(|v| {
            let (val, g) = ddg_points_as_tori::csg::ImplicitField::distance_gradient(&f, *v).unwrap();
            (val / (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt()).abs()
        })
        .collect();
    errs.sort_by(f64::total_cmp);
    assert!(errs[errs.len() / 2] < 0.05 * h, "median distance to F=0: {}", errs[errs.len() / 2]);
    // simplification keeps the volume and stays close
    let (sv, sf) = simplify_closed(&m.vertices, &m.faces, 1024);
    assert!(sf.len() <= 1024);
    let (r0, r1) = (check_mesh(&m.vertices, &m.faces).unwrap(), check_mesh(&sv, &sf).unwrap());
    assert!(r1.is_closed() && r1.degenerate_faces == 0);
    assert!((r1.signed_volume / r0.signed_volume - 1.0).abs() < 0.03, "{} vs {}", r1.signed_volume, r0.signed_volume);
    let (truth, _) = MeshTruth::new(&m.vertices, &m.faces, &TruthOptions::default()).unwrap();
    let worst = sv.iter().map(|v| truth.closest(*v).0).fold(0.0, f64::max);
    assert!(worst < 0.05, "simplified vertices up to {worst} from the original surface");
}

#[test]
fn blob_records_do_not_depend_on_the_thread_count() {
    let layout = Layout::default();
    let opts = BlobOptions::default();
    let mode = CloudMode::Scan(ScanOptions { spacing: (0.05, 0.1), ..ScanOptions::default() });
    let r1 = make_pool(Some(1)).unwrap().install(|| generate_blob_records(7, 0, 4, &mode, &layout, &opts).unwrap());
    let r8 = make_pool(Some(8)).unwrap().install(|| generate_blob_records(7, 0, 4, &mode, &layout, &opts).unwrap());
    assert_eq!(r1, r8);
}

#[test]
fn bvh_matches_brute_force_ray_casting() {
    let (v, f) = uv_sphere(24, 12, 0.8);
    let bvh = TriBvh::new(&v, &f).unwrap();
    let mut rng = Rng::new(4);
    for _ in 0..2000 {
        let o = [rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0)];
        let d = rng.unit_sphere();
        let got = bvh.intersect(o, d);
        // brute force Möller–Trumbore
        let mut best: Option<(f64, u32)> = None;
        for (fi, t) in f.iter().enumerate() {
            let (a, b, c) = (v[t[0] as usize], v[t[1] as usize], v[t[2] as usize]);
            let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let p = [d[1] * e2[2] - d[2] * e2[1], d[2] * e2[0] - d[0] * e2[2], d[0] * e2[1] - d[1] * e2[0]];
            let det = e1[0] * p[0] + e1[1] * p[1] + e1[2] * p[2];
            if det.abs() < 1e-14 {
                continue;
            }
            let s = [o[0] - a[0], o[1] - a[1], o[2] - a[2]];
            let u = (s[0] * p[0] + s[1] * p[1] + s[2] * p[2]) / det;
            let q = [s[1] * e1[2] - s[2] * e1[1], s[2] * e1[0] - s[0] * e1[2], s[0] * e1[1] - s[1] * e1[0]];
            let w = (d[0] * q[0] + d[1] * q[1] + d[2] * q[2]) / det;
            let t_ = (e2[0] * q[0] + e2[1] * q[1] + e2[2] * q[2]) / det;
            if u >= 0.0 && w >= 0.0 && u + w <= 1.0 && t_ > 0.0 && best.is_none_or(|b| t_ < b.0) {
                best = Some((t_, fi as u32));
            }
        }
        match (got, best) {
            (Some(g), Some(b)) => assert!((g.t - b.0).abs() < 1e-9),
            (None, None) => {}
            (g, b) => panic!("ray {o:?} {d:?}: bvh {g:?} vs brute {b:?}"),
        }
    }
}

#[test]
fn scanner_hits_lie_on_the_mesh_and_face_the_cameras() {
    let (m, _) = generate_blob(42, 3, &BlobOptions::default()).unwrap();
    let bvh = TriBvh::new(&m.vertices, &m.faces).unwrap();
    let (truth, _) = MeshTruth::new(&m.vertices, &m.faces, &TruthOptions::default()).unwrap();
    let opts = ScanOptions { spacing: (0.02, 0.06), ..ScanOptions::default() };
    for seed in 0..4 {
        let s = simulate_scan(&bvh, &opts, seed).unwrap();
        assert!(s.rays > 1000 && s.hits >= s.positions.len() && s.positions.len() >= 300 && s.positions.len() <= 2048, "{} rays {} hits {} kept", s.rays, s.hits, s.positions.len());
        assert_eq!(s.positions.len(), s.normals.len());
        for (p, n) in s.positions.iter().zip(&s.normals) {
            assert!(truth.closest(*p).0 < 1e-9, "hit off the mesh");
            assert!(((n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt() - 1.0).abs() < 1e-12);
        }
        // the camera sees the outside of a closed outward mesh: a point slightly along its
        // normal is outside (positive distance sign by the truth's winding)
        let outside = s.positions.iter().zip(&s.normals).step_by(11).filter(|(p, n)| truth.signed_distance([p[0] + 1e-4 * n[0], p[1] + 1e-4 * n[1], p[2] + 1e-4 * n[2]]).unwrap() > 0.0).count();
        let total = s.positions.iter().step_by(11).count();
        assert!(outside as f64 >= 0.97 * total as f64, "{outside}/{total} scan normals point outwards");
    }
    // determinism by seed
    let a = simulate_scan(&bvh, &opts, 5).unwrap();
    let b = simulate_scan(&bvh, &opts, 5).unwrap();
    let c = simulate_scan(&bvh, &opts, 6).unwrap();
    assert_eq!(a, b);
    assert_ne!(a, c);
    let a1 = make_pool(Some(1)).unwrap().install(|| simulate_scan(&bvh, &opts, 5).unwrap());
    let a8 = make_pool(Some(8)).unwrap().install(|| simulate_scan(&bvh, &opts, 5).unwrap());
    assert_eq!(a1, a8);
    assert_eq!(a1, a);
    // fibonacci directions are unit and spread
    let f = fibonacci_sphere(4);
    assert!(f.iter().all(|d| ((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - 1.0).abs() < 1e-12));
}

#[test]
fn farthest_point_sampling_is_greedy_and_thread_independent() {
    let mut rng = Rng::new(2);
    let pts: Vec<[f64; 3]> = (0..3000).map(|_| rng.unit_sphere()).collect();
    let fps = farthest_point_sampling(&pts, 100, 0);
    // brute-force greedy
    let mut chosen = vec![0usize];
    let mut mind = vec![f64::INFINITY; pts.len()];
    while chosen.len() < 100 {
        let q = pts[*chosen.last().unwrap()];
        let mut best = (-1.0, 0);
        for (i, p) in pts.iter().enumerate() {
            let d = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2);
            mind[i] = f64::min(mind[i], d);
            if mind[i] > best.0 {
                best = (mind[i], i);
            }
        }
        chosen.push(best.1);
        mind[best.1] = 0.0;
    }
    assert_eq!(fps, chosen);
    assert_eq!(fps, make_pool(Some(1)).unwrap().install(|| farthest_point_sampling(&pts, 100, 0)));
    // better separated than a random subset
    let sep = |ix: &[usize]| {
        let mut m = f64::INFINITY;
        for (a, &i) in ix.iter().enumerate() {
            for &j in &ix[a + 1..] {
                m = f64::min(m, (0..3).map(|k| (pts[i][k] - pts[j][k]).powi(2)).sum::<f64>());
            }
        }
        m.sqrt()
    };
    let random: Vec<usize> = (0..100).collect();
    assert!(sep(&fps) > 3.0 * sep(&random));
    assert_eq!(farthest_point_sampling(&pts, 5000, 0).len(), 3000);
    let p2: Vec<[f64; 2]> = pts.iter().map(|p| [p[0], p[1]]).collect();
    assert_eq!(farthest_point_sampling_2d(&p2, 50, 0).len(), 50);
}

#[test]
fn uniform_mesh_sampling_is_area_uniform() {
    let (v, f) = uv_sphere(64, 32, 1.0);
    let mut rng = Rng::new(3);
    let n = 40_000;
    let (p, nn) = sample_mesh_uniform(&v, &f, n, &mut rng).unwrap();
    // On a faceted sphere y is uniform up to the faceting; compare to the mesh area measure: the
    // fraction of points above the equator is 1/2 by symmetry, and each point lies on a face.
    let (truth, _) = MeshTruth::new(&v, &f, &TruthOptions::default()).unwrap();
    for (x, _) in p.iter().zip(&nn).step_by(50) {
        assert!(truth.closest(*x).0 < 1e-12);
    }
    let up = p.iter().filter(|x| x[1] > 0.0).count() as f64 / n as f64;
    assert!((up - 0.5).abs() < 4.5 * (0.25 / n as f64).sqrt());
    let ys: Vec<f64> = p.iter().map(|x| (x[1] + 1.0) / 2.0).collect();
    assert!(ks_uniform(ys) < 1.95 / (n as f64).sqrt() + 0.01); // facets perturb y slightly
    assert!(nn.iter().zip(&p).all(|(a, x)| a[0] * x[0] + a[1] * x[1] + a[2] * x[2] > 0.9));
    let c = center_and_scale(&v.iter().map(|x| [x[0] * 3.0 + 1.0, x[1] * 3.0, x[2] * 3.0]).collect::<Vec<_>>());
    assert!(c.iter().flat_map(|x| x.iter()).fold(0.0f64, |a, &b| a.max(b.abs())) <= 1.0 + 1e-12);
}

#[test]
fn training_npz_has_the_authors_layout() {
    let layout = Layout::default();
    let opts = BlobOptions::default();
    let recs = generate_blob_records(11, 0, 3, &CloudMode::Uniform, &layout, &opts).unwrap();
    let mut recs = recs;
    // a padded scanned shape
    let (m, _) = generate_blob(11, 100, &opts).unwrap();
    let (p, n) = sample_cloud(&m, &CloudMode::Scan(ScanOptions { spacing: (0.2, 0.3), ..ScanOptions::default() }), 2048, 1).unwrap();
    assert!(p.len() < 2048 && p.len() >= 20);
    recs.push(ShapeRecord { positions: p.clone(), normals: n, vertices: m.vertices.clone(), faces: m.faces.clone() });
    let dir = std::env::temp_dir().join(format!("pat-wp8-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("blobs_test.npz");
    write_training_npz(&path, &recs, &layout).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    let members = parse_npz_stored(&bytes).unwrap();
    let get = |name: &str| members.iter().find(|(n, _)| n == name).map(|(_, a)| a.clone()).unwrap();
    let (s, pp) = (4, 2048);
    assert_eq!(get("positions").shape, vec![s, pp, 3]);
    assert_eq!(get("positions").descr, "<f8");
    assert_eq!(get("normals").shape, vec![s, pp, 3]);
    assert_eq!(get("neighbors").shape, vec![s, pp, 64]);
    assert_eq!(get("neighbors").descr, "<i4");
    assert_eq!(get("masks").shape, vec![s, pp]);
    assert_eq!(get("masks").descr, "|b1");
    assert_eq!(get("mesh_vertices").shape, vec![s, 2048, 3]);
    assert_eq!(get("mesh_faces").shape, vec![s, 2048, 3]);
    assert_eq!(get("face_masks").shape, vec![s, 2048]);
    let masks = get("masks").data;
    let nb = get("neighbors").data;
    let pos = get("positions").data;
    let k = 64;
    for sh in 0..s {
        let nvalid = masks[sh * pp..(sh + 1) * pp].iter().filter(|&&b| b > 0.0).count();
        assert!(masks[sh * pp..sh * pp + nvalid].iter().all(|&b| b > 0.0), "mask is a prefix");
        for i in 0..nvalid {
            let row = &nb[(sh * pp + i) * k..(sh * pp + i + 1) * k];
            assert_eq!(row[0] as usize, i);
            let valid = row.iter().take_while(|&&j| j >= 0.0).count();
            assert_eq!(valid, k.min(nvalid));
            assert!(row[valid..].iter().all(|&j| j == -1.0));
            assert!(row[..valid].iter().all(|&j| (j as usize) < nvalid));
            // sorted by distance (after self)
            let d = |j: f64| (0..3).map(|c| (pos[(sh * pp + j as usize) * 3 + c] - pos[(sh * pp + i) * 3 + c]).powi(2)).sum::<f64>();
            assert!(row[1..valid].windows(2).all(|w| d(w[0]) <= d(w[1])));
        }
        let nf = recs[sh].faces.len();
        if nf < 2048 {
            assert_eq!(get("mesh_faces").data[(sh * 2048 + nf) * 3], -1.0);
        }
        assert!(get("mesh_faces").data[(sh * 2048 + nf - 1) * 3] >= 0.0);
    }
    // read back with the training reader when its feature is on
    #[cfg(feature = "train")]
    {
        let (shapes, rep) = ddg_points_as_tori::train::data::load_file(&path, &ddg_points_as_tori::train::data::DataOptions::default()).unwrap();
        assert_eq!(shapes.len(), 4);
        assert_eq!(rep.nonunit_normals, 0);
        assert_eq!(shapes[3].n, p.len());
        assert_eq!(shapes[0].n, 2048);
        assert_eq!(shapes[0].faces.len(), recs[0].faces.len());
    }
    std::fs::remove_dir_all(&dir).ok();
    // invalid shapes are errors, not panics
    let bad = ShapeRecord { positions: vec![[0.0; 3]; 3000], normals: vec![[0.0, 0.0, 1.0]; 3000], vertices: vec![], faces: vec![] };
    assert!(write_training_npz(&std::env::temp_dir().join("pat-wp8-bad.npz"), &[bad], &layout).is_err());
}

#[test]
fn marching_cubes_sphere_mesh_for_simplification_and_scan() {
    // scanning a sphere mesh: all hit radii equal the sphere radius within faceting
    let (v, f) = uv_sphere(48, 24, 0.9);
    let bvh = TriBvh::new(&v, &f).unwrap();
    let s = simulate_scan(&bvh, &ScanOptions { spacing: (0.05, 0.06), ..ScanOptions::default() }, 3).unwrap();
    assert!(s.positions.iter().all(|p| {
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        r > 0.9 * (PI / 24.0).cos() - 1e-9 && r < 0.9 + 1e-9
    }));
    assert!(s.normals.iter().zip(&s.positions).all(|(n, p)| n[0] * p[0] + n[1] * p[1] + n[2] * p[2] > 0.0));
    // MC of a sphere field then simplification to 200 faces keeps a closed manifold
    struct S;
    impl ddg_points_as_tori::csg::ImplicitField for S {
        fn distance(&self, x: [f64; 3]) -> Result<f64, ddg_points_as_tori::PatError> {
            Ok((x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt() - 0.7)
        }
        fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), ddg_points_as_tori::PatError> {
            let l = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
            Ok((l - 0.7, [x[0] / l, x[1] / l, x[2] / l]))
        }
        fn bbox(&self) -> [[f64; 3]; 2] {
            [[-1.0; 3], [1.0; 3]]
        }
    }
    let (m, _) = marching_cubes(&S, &McOptions { resolution: Resolution::Longest(32), band: Band::Dense, ..McOptions::default() }).unwrap();
    let (sv, sf) = simplify_closed(&m.positions, &m.faces, 200);
    assert!(sf.len() <= 200);
    let r = check_mesh(&sv, &sf).unwrap();
    assert!(r.is_closed() && (r.signed_volume / (4.0 / 3.0 * PI * 0.343) - 1.0).abs() < 0.1, "{r:?}");
}

#[test]
fn two_d_generators() {
    // polygon distance against brute force; sign by area
    let mut rng = Rng::new(6);
    let poly = random_polygon(&mut rng, 10).unwrap();
    let cloud = poly.sample_uniform(500, &mut rng);
    for (p, n) in cloud.positions.iter().zip(&cloud.normals) {
        let (d, g) = poly.signed_distance(*p);
        assert!(d.abs() < 1e-12);
        assert!((g[0] * n[0] + g[1] * n[1] - 1.0).abs() < 1e-9 || d.abs() < 1e-12);
        // outward normals: a step outwards is outside
        assert!(poly.signed_distance([p[0] + 1e-6 * n[0], p[1] + 1e-6 * n[1]]).0 > 0.0);
        assert!(poly.signed_distance([p[0] - 1e-6 * n[0], p[1] - 1e-6 * n[1]]).0 < 0.0);
    }
    // noise blob
    let params = Blob2dParams { alpha: 0.4, scale: 1.2, octaves: 3, noise_base: 5 };
    let blob = noise_blob(&params, 256).unwrap();
    let bb = blob.bbox();
    assert!(bb.iter().flatten().fold(0.0f64, |a, &c| a.max(c.abs())) <= 1.0 + 1e-9);
    assert!(blob.signed_distance([0.0, 0.0]).0 < 0.0 || blob.signed_distance([0.1, 0.0]).0 < 0.0);
    // scanner: hits on the curve, normals facing, determinism
    let opts = Scan2dOptions::default();
    let a = simulate_scan_2d(&blob, &opts, 1).unwrap();
    let b = simulate_scan_2d(&blob, &opts, 1).unwrap();
    assert_eq!(a, b);
    assert!(a.positions.len() > 50 && a.positions.len() <= 512, "{}", a.positions.len());
    for (p, n) in a.positions.iter().zip(&a.normals) {
        assert!(blob.signed_distance(*p).0.abs() < 1e-9);
        assert!(((n[0] * n[0] + n[1] * n[1]).sqrt() - 1.0).abs() < 1e-12);
    }
    // analytic instance: boundary samples have exact normals and zero distance
    let inst = Instance2d::Analytic(Shape2d::new(Kind::Moon, &[0.5, 0.6, 0.3]).unwrap());
    let c = inst.sample_uniform(300, &mut rng).unwrap();
    for (p, n) in c.positions.iter().zip(&c.normals) {
        let (d, g) = inst.signed_distance(*p);
        assert!(d.abs() < 1e-11 && (g[0] * n[0] + g[1] * n[1] - 1.0).abs() < 1e-12);
    }
    let sc = inst.scan(&opts, 2).unwrap();
    assert!(sc.positions.iter().all(|p| inst.signed_distance(*p).0.abs() < 5e-3));
    // random instances are deterministic and mixed
    let mut kinds = [0usize; 3];
    for i in 0..60 {
        let x = random_instance_2d(1, i).unwrap();
        let y = random_instance_2d(1, i).unwrap();
        assert_eq!(x.signed_distance([0.2, 0.1]), y.signed_distance([0.2, 0.1]));
        kinds[match x {
            Instance2d::Analytic(_) => 0,
            Instance2d::Polygon(_) => 1,
            Instance2d::Blob(..) => 2,
        }] += 1;
    }
    assert!(kinds.iter().all(|&k| k >= 8), "{kinds:?}");
    assert!(Curve2d::polygon(&[[0.0, 0.0], [1.0, 0.0]]).is_err());
}
