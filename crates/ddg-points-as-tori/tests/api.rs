//! API items of SPEC.md §18 (WP12/WP13 resolutions): root re-exports, public `Morph::levels`,
//! `morphology` returning its cloud, field accessors for the cloud and its hash, the sign-rule
//! marker of the tori file, `smooth = 0`, the evaluation report of marching cubes, and the
//! progress/cancel hook.

mod common;

use common::{bits, ellipsoid_cloud, lcg};
use ddg_points_as_tori::apps::{Morph, MorphOptions, SampleOptions, morphology};
use ddg_points_as_tori::csg::{Csg, ImplicitField};
use ddg_points_as_tori::mc::{McOptions, Resolution, marching_cubes};
use ddg_points_as_tori::{BuildOptions, EigenRule, EvalOut, GpuMatmul, Model, PatError, PointsAsTori, Progress, SignRule, Tori, TorusField};
use std::sync::atomic::{AtomicBool, Ordering};

fn built(n: usize, rule: SignRule) -> (ddg_points_as_tori::PointCloud, TorusField) {
    let cloud = ellipsoid_cloud(n, 3, [0.6, 0.45, 0.3]);
    let mut o = BuildOptions::default().with_threads(Some(8));
    o.fit.sign_rule = rule;
    let (f, _) = PointsAsTori::build(&cloud, &Model::reference().unwrap(), &o).unwrap();
    (cloud, f)
}

#[test]
fn root_reexports() {
    // compile-time: the types are at the crate root
    let _ = (EigenRule::Robust, GpuMatmul::Split3, SignRule::Geometric);
    assert_eq!(Morph::Close(0.1).levels().unwrap(), vec![0.1, -0.1]);
    assert_eq!(Morph::Open(0.1).levels().unwrap(), vec![-0.1, 0.1]);
    assert_eq!(Morph::Erode(0.2).levels().unwrap(), vec![-0.2]);
    assert!(Morph::Dilate(-1.0).levels().is_err());
    assert!(Morph::Dilate(f64::NAN).levels().is_err());
}

#[test]
fn field_keeps_its_cloud_and_hash() {
    let (cloud, f) = built(500, SignRule::Reference);
    assert_eq!(f.cloud(), &cloud);
    assert_eq!(f.cloud_hash(), cloud.hash());
    assert_eq!(f.sign_rule(), Some(SignRule::Reference));
    let (_, g) = built(500, SignRule::Geometric);
    assert_eq!(g.sign_rule(), Some(SignRule::Geometric));
}

#[test]
fn tori_file_marks_the_sign_rule_and_reads_old_files() {
    let (cloud, f) = built(500, SignRule::Geometric);
    // new file
    let bytes = f.to_bytes();
    let (t, h, rule) = Tori::from_bytes_with(&bytes, Some(cloud.hash())).unwrap();
    assert_eq!((h, rule), (cloud.hash(), Some(SignRule::Geometric)));
    assert_eq!(&t, f.tori());
    // the old entry points still read it, and still write "not recorded"
    assert_eq!(Tori::from_bytes(&bytes, None).unwrap().0, t);
    let old = f.tori().to_bytes(cloud.hash());
    assert_eq!(&old[12..16], &[0, 0, 0, 0]);
    assert_eq!(Tori::from_bytes_with(&old, None).unwrap().2, None);
    // identical but for the marker
    assert_eq!(old[..12], bytes[..12]);
    assert_eq!(old[16..], bytes[16..]);
    // an unknown marker is an error, a wrong cloud too
    let mut bad = bytes.clone();
    bad[12] = 9;
    assert!(matches!(Tori::from_bytes_with(&bad, None), Err(PatError::Format { .. })));
    assert!(Tori::from_bytes_with(&bytes, Some(cloud.hash() ^ 1)).is_err());
    // Reference marker
    let r = f.tori().to_bytes_with(cloud.hash(), Some(SignRule::Reference));
    assert_eq!(Tori::from_bytes_with(&r, None).unwrap().2, Some(SignRule::Reference));
}

#[test]
fn morphology_returns_the_sampled_cloud() {
    let (_, f) = built(800, SignRule::Reference);
    let opts = MorphOptions {
        sample: SampleOptions { count: Some(500), mc: McOptions { resolution: Resolution::Longest(24), threads: Some(4), ..Default::default() }, threads: Some(4), ..Default::default() },
        build: BuildOptions::default().with_threads(Some(4)),
        paper_eval: true,
    };
    let (g, cloud, rep) = morphology(&f, &Model::reference().unwrap(), Morph::Dilate(0.03), &opts).unwrap();
    assert_eq!(g.cloud(), &cloud);
    assert_eq!(cloud.len(), 500);
    assert_eq!(rep.steps.len(), 1);
}

struct Ball(f64);
impl ImplicitField for Ball {
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

#[test]
fn smooth_zero_is_the_sharp_boolean() {
    let a = || Box::new(Csg::leaf(Ball(0.5)));
    let b = || Box::new(Csg::leaf(Ball(0.4)).offset(-0.1));
    let mut s = 3;
    let xs: Vec<[f64; 3]> = (0..200).map(|_| [2.0 * lcg(&mut s) - 1.0, 2.0 * lcg(&mut s) - 1.0, 2.0 * lcg(&mut s) - 1.0]).collect();
    for (smooth, sharp) in [
        (Csg::SmoothUnion(a(), b(), 0.0), Csg::Union(a(), b())),
        (Csg::SmoothIntersect(a(), b(), 0.0), Csg::Intersect(a(), b())),
        (Csg::SmoothSubtract(a(), b(), 0.0), Csg::Subtract(a(), b())),
    ] {
        for x in &xs {
            assert_eq!(smooth.distance_gradient(*x).unwrap(), sharp.distance_gradient(*x).unwrap());
        }
    }
    assert!(Csg::SmoothUnion(a(), b(), -1e-3).distance([0.0; 3]).is_err());
}

#[test]
fn marching_cubes_reports_the_evaluations_of_the_fitted_fields() {
    let (_, f) = built(600, SignRule::Reference);
    let o = McOptions { resolution: Resolution::Longest(20), threads: Some(4), ..Default::default() };
    let (_, rep) = marching_cubes(&f, &o).unwrap();
    let e = rep.eval.expect("a TorusField reports its evaluations");
    assert!(e.queries >= rep.grid_values && e.terms >= 32 * e.queries as u64 - 1, "{e:?} {}", rep.grid_values);
    // a second run reports its own, not the sum
    let (_, rep2) = marching_cubes(&f, &o).unwrap();
    assert_eq!(rep2.eval.unwrap().queries, e.queries);
    assert_eq!(f.totals().queries, 2 * e.queries);
    // a tree sums its leaves; analytic fields report nothing
    let tree = Csg::field(built(600, SignRule::Reference).1).union(Csg::leaf(Ball(0.2)));
    assert!(marching_cubes(&tree, &o).unwrap().1.eval.is_some());
    assert!(marching_cubes(&Ball(0.5), &o).unwrap().1.eval.is_none());
}

#[test]
fn progress_counts_and_changes_nothing() {
    let cloud = ellipsoid_cloud(1500, 5, [0.6, 0.45, 0.3]);
    let model = Model::reference().unwrap();
    let plain = PointsAsTori::build(&cloud, &model, &BuildOptions::default().with_threads(Some(8))).unwrap().0;
    let p = Progress::new();
    assert_eq!((p.stage(), p.done(), p.total(), p.is_cancelled()), ("idle", 0, 0, false));
    let mut o = BuildOptions::default().with_threads(Some(8));
    o.progress = Some(p.clone());
    o.eval.progress = Some(p.clone());
    assert_eq!(o, o.clone());
    let f = PointsAsTori::build(&cloud, &model, &o).unwrap().0;
    assert_eq!(p.stage(), "tree");
    assert_eq!((p.done(), p.total()), (1, 1));
    assert_eq!(f.tori(), plain.tori());
    // evaluation: counts the queries, same bits
    let qs: Vec<[f64; 3]> = (0..2000).map(|i| [0.01 * (i % 90) as f64 - 0.4, 0.01 * (i % 70) as f64 - 0.3, 0.001 * i as f64 - 1.0]).collect();
    let (mut a, mut b) = (vec![0.0; qs.len()], vec![0.0; qs.len()]);
    f.evaluate(&qs, EvalOut { distance: Some(&mut a), gradient: None, laplacian: None }).unwrap();
    assert_eq!((p.stage(), p.done(), p.total(), p.fraction()), ("evaluate", 2000, 2000, 1.0));
    plain.evaluate(&qs, EvalOut { distance: Some(&mut b), gradient: None, laplacian: None }).unwrap();
    assert_eq!(bits(&a), bits(&b));
}

#[test]
fn cancel_stops_builds_and_evaluations() {
    let cloud = ellipsoid_cloud(6000, 5, [0.6, 0.45, 0.3]);
    let model = Model::reference().unwrap();
    // before the call
    let p = Progress::new();
    p.cancel();
    let mut o = BuildOptions::default().with_threads(Some(8));
    o.progress = Some(p.clone());
    assert_eq!(PointsAsTori::build(&cloud, &model, &o).unwrap_err(), PatError::Cancelled);
    // during the network stage, from another thread
    p.reset();
    assert!(!p.is_cancelled());
    let done = AtomicBool::new(false);
    let res = std::thread::scope(|s| {
        s.spawn(|| {
            while !done.load(Ordering::Relaxed) {
                if p.stage() == "network" && p.done() > 0 {
                    p.cancel();
                    return;
                }
                std::thread::yield_now();
            }
        });
        let r = PointsAsTori::build(&cloud, &model, &o);
        done.store(true, Ordering::Relaxed);
        r
    });
    assert_eq!(res.unwrap_err(), PatError::Cancelled);
    assert!(p.done() < p.total(), "stopped early: {} of {}", p.done(), p.total());
    // an evaluation
    let (_, f) = built(500, SignRule::Reference);
    let q = Progress::new();
    let g = TorusField::new(f.cloud(), f.tori().clone(), ddg_points_as_tori::EvalOptions { progress: Some(q.clone()), threads: Some(4), ..Default::default() }).unwrap();
    q.cancel();
    let mut d = vec![0.0; 10];
    assert_eq!(g.evaluate(&[[0.0; 3]; 10], EvalOut { distance: Some(&mut d), gradient: None, laplacian: None }).unwrap_err(), PatError::Cancelled);
    q.reset();
    g.evaluate(&[[0.1; 3]; 10], EvalOut { distance: Some(&mut d), gradient: None, laplacian: None }).unwrap();
}
