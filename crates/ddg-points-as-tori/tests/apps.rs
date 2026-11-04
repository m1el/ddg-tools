//! Level-set sampling, offsets of fitted fields and morphology (SPEC.md §8.1–8.3). Run with `--release`.

mod common;

use ddg_points_as_tori::apps::{Morph, MorphOptions, SampleOptions, morphology, sample_level_set};
use ddg_points_as_tori::csg::{Csg, ImplicitField};
use ddg_points_as_tori::mc::{McOptions, Resolution, marching_cubes};
use ddg_points_as_tori::{BuildOptions, EvalOptions, Model, NormalPolicy, PatError, PointCloud, PointsAsTori, TorusField};

struct Sphere(f64);
impl ImplicitField for Sphere {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        Ok((x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt() - self.0)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let l = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        Ok((l - self.0, [x[0] / l, x[1] / l, x[2] / l]))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        [[-1.0; 3], [1.0; 3]]
    }
}

fn fibonacci_sphere(n: usize, r: f64) -> PointCloud {
    let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    let (mut p, mut nr) = (Vec::new(), Vec::new());
    for i in 0..n {
        let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
        let rr = (1.0 - z * z).sqrt();
        let t = golden * i as f64;
        let u = [rr * t.cos(), rr * t.sin(), z];
        p.push(u.map(|c| c * r));
        nr.push(u);
    }
    PointCloud::new(p, nr, NormalPolicy::RequireUnit { tol: 1e-9 }).unwrap().0
}

fn build(cloud: &PointCloud, threads: usize) -> TorusField {
    let model = Model::reference().unwrap();
    let mut b = BuildOptions::default().with_threads(Some(threads));
    b.eval = TorusField::paper_options(cloud, Some(threads)).unwrap().0;
    PointsAsTori::build(cloud, &model, &b).unwrap().0
}

#[test]
fn level_set_samples_have_small_residuals_and_unit_normals() {
    let s = Sphere(0.5);
    for level in [-0.2, 0.0, 0.15] {
        let opts = SampleOptions {
            count: Some(2000),
            mc: McOptions { resolution: Resolution::Longest(32), threads: Some(4), ..Default::default() },
            threads: Some(4),
            ..Default::default()
        };
        let (p, n, rep) = sample_level_set(&s, level, &opts).unwrap();
        assert_eq!((p.len(), rep.produced, rep.rejected), (2000, 2000, 0), "level {level}");
        assert!(rep.residual_max < 1e-12, "{}", rep.residual_max);
        assert!(rep.residual_before_max > 1e-6, "MC alone should not be this exact");
        for (x, nn) in p.iter().zip(&n) {
            assert!((s.distance(*x).unwrap() - level).abs() < 1e-12);
            assert!((nn.iter().map(|c| c * c).sum::<f64>() - 1.0).abs() < 1e-12);
            let r = x.iter().map(|c| c * c).sum::<f64>().sqrt();
            assert!((0..3).all(|i| (nn[i] - x[i] / r).abs() < 1e-12));
        }
        // area-uniform: the octant counts of a sphere are equal up to sampling noise
        let mut oct = [0usize; 8];
        for x in &p {
            oct[(x[0] > 0.0) as usize + 2 * (x[1] > 0.0) as usize + 4 * (x[2] > 0.0) as usize] += 1;
        }
        assert!(oct.iter().all(|c| (*c as f64 - 250.0).abs() < 60.0), "{oct:?}");
        // deterministic and thread independent
        let (p2, _, _) = sample_level_set(&s, level, &SampleOptions { threads: Some(1), mc: McOptions { threads: Some(1), ..opts.mc.clone() }, ..opts.clone() }).unwrap();
        assert_eq!(p, p2);
    }
    // the empty level set is an error naming the level
    let e = sample_level_set(&s, 10.0, &SampleOptions { count: Some(10), ..Default::default() }).unwrap_err();
    assert!(e.to_string().contains("empty"), "{e}");
    assert!(sample_level_set(&s, f64::NAN, &SampleOptions::default()).is_err());
    assert!(sample_level_set(&s, 0.0, &SampleOptions { newton_iters: 0, ..Default::default() }).is_err());
}

#[test]
fn torus_field_marching_cubes_and_offset() {
    let cloud = fibonacci_sphere(512, 0.5);
    let f = build(&cloud, 8);
    let (m, rep) = marching_cubes(&f, &McOptions { resolution: Resolution::Longest(48), threads: Some(8), ..Default::default() }).unwrap();
    assert!(!rep.bound_violated, "max slope {}", rep.max_edge_slope);
    assert!(rep.max_edge_slope < 2.0, "{}", rep.max_edge_slope);
    match rep.mesh {
        ddg_points_as_tori::mc::MeshCheck::Manifold { boundary_edges, euler_characteristic, .. } => assert_eq!((boundary_edges, euler_characteristic), (0, 2)),
        ref o => panic!("{o:?}"),
    }
    // bitwise independent of the marching cubes thread count (the field's own pool is fixed)
    let (m1, _) = marching_cubes(&f, &McOptions { resolution: Resolution::Longest(48), threads: Some(1), ..Default::default() }).unwrap();
    assert_eq!(m1.faces, m.faces);
    assert!(m1.positions.iter().zip(&m.positions).all(|(a, b)| a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())));
    let r: Vec<f64> = m.positions.iter().map(|p| p.iter().map(|c| c * c).sum::<f64>().sqrt()).collect();
    let (lo, hi) = (r.iter().cloned().fold(f64::INFINITY, f64::min), r.iter().cloned().fold(0.0, f64::max));
    assert!(lo > 0.48 && hi < 0.52, "radius range {lo}..{hi}");
    // offset by c: the level set phi = c is a sphere of radius about 0.5 + c (SPEC 8.1, no refit)
    let (m2, _) = marching_cubes(&f, &McOptions { isovalue: 0.1, resolution: Resolution::Longest(48), threads: Some(8), ..Default::default() }).unwrap();
    let r2: Vec<f64> = m2.positions.iter().map(|p| p.iter().map(|c| c * c).sum::<f64>().sqrt()).collect();
    let mean = r2.iter().sum::<f64>() / r2.len() as f64;
    assert!((mean - 0.6).abs() < 0.02, "{mean}");
    // a CSG of fitted fields: gradient of the union is the gradient of the active argument
    let other = Csg::leaf(Sphere(0.3));
    let u = Csg::field(f).union(other);
    let x = [0.45, 0.1, 0.0];
    let (d, g) = u.distance_gradient(x).unwrap();
    assert!(d.abs() < 0.1);
    assert!((g.iter().map(|c| c * c).sum::<f64>().sqrt() - 1.0).abs() < 0.2);
}

fn hausdorff_to_field(points: &[[f64; 3]], f: &dyn ImplicitField) -> f64 {
    points.iter().map(|p| f.distance(*p).unwrap().abs()).fold(0.0, f64::max)
}

#[test]
fn morphology_round_trip_on_a_sphere() {
    let cloud = fibonacci_sphere(512, 0.5);
    let f0 = build(&cloud, 8);
    let model = Model::reference().unwrap();
    let opts = MorphOptions {
        build: BuildOptions::default().with_threads(Some(8)),
        sample: SampleOptions {
            threads: Some(8),
            mc: McOptions { resolution: Resolution::Longest(48), threads: Some(8), ..Default::default() },
            ..Default::default()
        },
        ..Default::default()
    };
    let c = 0.05;
    // erosion / dilation shift the radius by c
    let (er, _cl, rep) = morphology(&f0, &model, Morph::Erode(c), &opts).unwrap();
    assert_eq!(rep.steps.len(), 1);
    assert!(rep.steps[0].sample.produced >= 500, "{}", rep.steps[0].sample.produced);
    let probe = fibonacci_sphere(400, 0.5 - c);
    assert!(hausdorff_to_field(probe.positions(), &er) < 0.02, "erode: {}", hausdorff_to_field(probe.positions(), &er));
    let (di, _, _) = morphology(&f0, &model, Morph::Dilate(c), &opts).unwrap();
    let probe = fibonacci_sphere(400, 0.5 + c);
    assert!(hausdorff_to_field(probe.positions(), &di) < 0.02, "dilate: {}", hausdorff_to_field(probe.positions(), &di));
    // opening and closing of a convex shape roughly preserve it
    for op in [Morph::Open(c), Morph::Close(c)] {
        let (g, _cl, rep) = morphology(&f0, &model, op, &opts).unwrap();
        assert_eq!(rep.steps.len(), 2);
        let err = hausdorff_to_field(cloud.positions(), &g);
        assert!(err < 0.03, "{op:?}: original points are {err} from the result's surface");
        let (m, _) = marching_cubes(&g, &McOptions { resolution: Resolution::Longest(48), threads: Some(8), ..Default::default() }).unwrap();
        let exact = 4.0 / 3.0 * std::f64::consts::PI * 0.5f64.powi(3);
        assert!((m.signed_volume() - exact).abs() < 0.15 * exact, "{op:?}: volume {} vs {exact}", m.signed_volume());
    }
    assert!(morphology(&f0, &model, Morph::Open(-1.0), &opts).is_err());
    assert!(morphology(&f0, &model, Morph::Erode(f64::NAN), &opts).is_err());
    // eroding past the radius leaves an empty level set: an error, not a panic
    assert!(morphology(&f0, &model, Morph::Erode(2.0), &opts).is_err());
}

#[test]
fn morphology_round_trip_on_the_golden_torus() {
    let Some(root) = common::golden_root() else { return };
    let g = common::Golden::open(&root, "synth_torus");
    let f0 = build(&g.cloud, 8);
    let model = Model::reference().unwrap();
    let opts = MorphOptions {
        build: BuildOptions::default().with_threads(Some(8)),
        sample: SampleOptions {
            threads: Some(8),
            mc: McOptions { resolution: Resolution::Longest(64), threads: Some(8), ..Default::default() },
            ..Default::default()
        },
        ..Default::default()
    };
    // the tube radius is 0.15; closing by 0.03 keeps a convex-ish torus nearly unchanged
    let (c, _, _) = morphology(&f0, &model, Morph::Close(0.03), &opts).unwrap();
    let err = hausdorff_to_field(g.cloud.positions(), &c);
    assert!(err < 0.04, "closed torus: {err}");
    let (o, _, _) = morphology(&f0, &model, Morph::Open(0.03), &opts).unwrap();
    let err = hausdorff_to_field(g.cloud.positions(), &o);
    assert!(err < 0.04, "opened torus: {err}");
    let _ = EvalOptions::default();
}
