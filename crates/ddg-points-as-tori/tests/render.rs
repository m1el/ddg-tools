//! WP11: the sphere tracer (SPEC.md §8.5). Analytic sphere, CSG, a fitted field; determinism at
//! several thread counts; bad input returns errors.

mod common;

use common::ellipsoid_cloud;
use ddg_points_as_tori::PatError;
use ddg_points_as_tori::csg::{Csg, ImplicitField};
use ddg_points_as_tori::render::{Camera, CutPlane, Environment, Epsilon, Tracer, TraceOptions, status, trace_ray, RayJob, TraceParams};
use ddg_points_as_tori::{BuildOptions, Model, PointsAsTori};

struct Sphere {
    c: [f64; 3],
    r: f64,
}

impl ImplicitField for Sphere {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        Ok(self.distance_gradient(x)?.0)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let d = [x[0] - self.c[0], x[1] - self.c[1], x[2] - self.c[2]];
        let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        Ok((l - self.r, if l > 0.0 { [d[0] / l, d[1] / l, d[2] / l] } else { [0.0, 1.0, 0.0] }))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        let p = self.r * 1.2;
        [[self.c[0] - p, self.c[1] - p, self.c[2] - p], [self.c[0] + p, self.c[1] + p, self.c[2] + p]]
    }
}

fn cam(w: usize, h: usize) -> Camera {
    Camera::look_at([0.0, 0.5, 5.0], [0.0; 3], [0.0, 1.0, 0.0], 30.0, w, h).unwrap()
}

#[test]
fn sphere_hit_position_and_report() {
    let s = Sphere { c: [0.0; 3], r: 1.0 };
    let job = RayJob { ro: [0.0, 0.0, 5.0], rd: [0.0, 0.0, -1.0], t0: 0.0, t1: 10.0 };
    let p = TraceParams { level: 0.0, safety: 0.9, g_min: 0.1, eps: 1e-6, max_steps: 128 };
    let r = trace_ray(&s, &job, &p);
    assert_eq!(r.status, status::HIT);
    assert!((r.t - 4.0).abs() < 1e-5, "{}", r.t);
    // an offset level 0.5 hits at radius 1.5
    let r = trace_ray(&s, &job, &TraceParams { level: 0.5, ..p });
    assert!((r.t - 3.5).abs() < 1e-5);
    // a ray that misses
    let off = RayJob { ro: [3.0, 0.0, 5.0], ..job };
    assert_eq!(trace_ray(&s, &off, &p).status, status::MISS);
    // step limit
    let r = trace_ray(&s, &job, &TraceParams { max_steps: 2, ..p });
    assert_eq!(r.status, status::STEP_LIMIT);

    let tr = Tracer::new(TraceOptions { threads: Some(2), ..Default::default() }).unwrap();
    let (img, rep) = tr.render(&s, &cam(64, 48)).unwrap();
    assert_eq!(img.rgb.len(), 64 * 48 * 3);
    assert_eq!(rep.rays, 64 * 48);
    assert!(rep.hits > 300 && rep.hits < 3000, "{rep:?}");
    assert_eq!(rep.step_limit, 0);
    assert!(rep.steps_max <= 128 && rep.steps_mean >= 1.0);
    // the centre pixel is the sphere (not background)
    let c = (24 * 64 + 32) * 3;
    assert_ne!(&img.rgb[c..c + 3], &[255, 255, 255]);
    assert_eq!(img.alpha[24 * 64 + 32], 255);
    // PNG round trip is a valid PNG
    let png = img.to_png_bytes(true).unwrap();
    assert_eq!(&png[1..4], b"PNG");
}

#[test]
fn csg_and_cut_plane_render() {
    let a = Csg::leaf(Sphere { c: [-0.4, 0.0, 0.0], r: 0.8 });
    let b = Csg::leaf(Sphere { c: [0.4, 0.0, 0.0], r: 0.8 });
    let scene = a.union(b).offset(0.1);
    let opts = TraceOptions { cut_plane: Some(CutPlane::default()), threads: Some(2), ..Default::default() };
    let tr = Tracer::new(opts).unwrap();
    let (img, rep) = tr.render(&scene, &cam(80, 60)).unwrap();
    assert!(rep.hits > 100);
    assert!(rep.cut_plane_pixels > 0, "{rep:?}");
    assert!(img.alpha.iter().any(|a| *a == 255));
}

#[test]
fn constant_sky_and_demo_options() {
    let s = Sphere { c: [0.0; 3], r: 1.0 };
    for opts in [
        TraceOptions { environment: Environment::Constant([0.5, 0.6, 0.7]), ..Default::default() },
        TraceOptions { threads: Some(1), ..TraceOptions::demo() },
    ] {
        let tr = Tracer::new(opts).unwrap();
        let cam = Camera::demo([0.0; 3], 0.5, 40, 30).unwrap();
        let (_, rep) = tr.render(&s, &cam).unwrap();
        assert!(rep.hits > 0, "{rep:?}");
    }
}

#[test]
fn bad_input_is_an_error() {
    let s = Sphere { c: [0.0; 3], r: 1.0 };
    assert!(Camera::look_at([0.0; 3], [0.0; 3], [0.0, 1.0, 0.0], 30.0, 8, 8).is_err());
    assert!(Camera::look_at([0.0, 0.0, 1.0], [0.0; 3], [0.0, 0.0, 1.0], 30.0, 8, 8).is_err());
    assert!(Camera::look_at([0.0, 0.0, 1.0], [0.0; 3], [0.0, 1.0, 0.0], 0.0, 8, 8).is_err());
    assert!(Camera::look_at([0.0, 0.0, 1.0], [0.0; 3], [0.0, 1.0, 0.0], 30.0, 0, 8).is_err());
    assert!(Camera::look_at([0.0, f64::NAN, 1.0], [0.0; 3], [0.0, 1.0, 0.0], 30.0, 8, 8).is_err());
    assert!(Camera::look_at([0.0, 0.0, 1.0], [0.0; 3], [0.0, 1.0, 0.0], 30.0, 1 << 20, 1 << 20).is_err());
    assert!(Tracer::new(TraceOptions { safety: 0.0, ..Default::default() }).is_err());
    assert!(Tracer::new(TraceOptions { safety: 1.5, ..Default::default() }).is_err());
    assert!(Tracer::new(TraceOptions { max_steps: 0, ..Default::default() }).is_err());
    assert!(Tracer::new(TraceOptions { eps_hit: Epsilon::Absolute(f64::NAN), ..Default::default() }).is_err());
    assert!(Tracer::new(TraceOptions { level: f64::INFINITY, ..Default::default() }).is_err());
    assert!(Tracer::new(TraceOptions { threads: Some(0), ..Default::default() }).unwrap().render(&s, &cam(8, 8)).is_err());
    assert!(Tracer::new(TraceOptions { environment: Environment::Constant([-1.0, 0.0, 0.0]), ..Default::default() }).is_err());
    assert!(Tracer::new(TraceOptions { cut_plane: Some(CutPlane { axis: 3, ..Default::default() }), ..Default::default() }).is_err());
    assert!(Tracer::new(TraceOptions { device: ddg_points_as_tori::Device::Cuda { ordinal: 0 }, ..Default::default() }).unwrap().render(&s, &cam(8, 8)).is_err());
    // a field with an empty box
    let flat = Sphere { c: [0.0; 3], r: 0.0 };
    assert!(Tracer::new(TraceOptions::default()).unwrap().render(&flat, &cam(8, 8)).is_err());
}

/// Bitwise identical image and report counters at 1 (twice), 3 and 8 threads, for a fitted field
/// with an offset and a Boolean, a cut plane and a shadow.
#[test]
fn determinism_across_threads() {
    let model = Model::reference().unwrap();
    let build = |seed, axes| {
        let cloud = ellipsoid_cloud(300, seed, axes);
        let (f, _) = PointsAsTori::build(&cloud, &model, &BuildOptions::default().with_threads(Some(4))).unwrap();
        f
    };
    let a = build(3, [1.0, 0.7, 0.5]);
    let b = build(4, [0.5, 0.9, 0.6]);
    let scene = Csg::field(a).offset(0.05).subtract(Csg::field(b));
    let cam = cam(48, 36);
    let run = |threads: usize| {
        let tr = Tracer::new(TraceOptions { threads: Some(threads), cut_plane: Some(CutPlane { axis: 0, ..Default::default() }), ..Default::default() }).unwrap();
        let (img, rep) = tr.render(&scene, &cam).unwrap();
        (img, (rep.hits, rep.misses, rep.clipped, rep.steps_total, rep.steps_max, rep.shadow_rays, rep.shadow_hits, rep.cut_plane_pixels, rep.step_limit))
    };
    let base = run(1);
    assert!(base.1.0 > 20, "no hits: {:?}", base.1);
    assert!(base == run(1), "two runs differ");
    for t in [3, 8] {
        assert!(base == run(t), "{t} threads differ");
    }
}
